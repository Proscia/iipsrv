// Member functions for CZIImage.h

/*  IIP Server: Tiled Pyramidal CZI handler

    Copyright (C) 2000-2019 Ruven Pillay.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

*/


#include "CZIImage.h"
#include <sstream>
#include <iostream>
#include <fstream>

#include "CziUtils.h"

extern std::ofstream logfile;

using namespace std;

void CZIImage::openImage()
{
  // Insist that the czi_reader and tile_buf be NULL
  if( czi_reader || tile_buf ){
    throw file_error( "CZIImage::openImage(): czi_reader or tile_buf is not NULL" );
  }

  string filename = getFileName( currentX, currentY );

  // Update our timestamp
  updateTimestamp( filename );

  // Try to allocate and open a CZI-reader.
  // See:  CZICmd/execute.cpp -- class CExecuteBase -- CreateAndOpenCziReader().
  std::wstring wide_filename = std::wstring(filename.begin(), filename.end());
  auto stream = libCZI::CreateStreamFromFile(wide_filename.c_str());
  czi_reader = libCZI::CreateCZIReader();
  czi_reader->Open(stream);

  // Load our metadata if not already loaded
  if( bpc == 0 ) loadImageInfo( currentX, currentY );

  // Insist on a tiled image
  if( (tile_width == 0) && (tile_height == 0) ){
    throw file_error( "CZI image is not tiled" );
  }

  isSet = true;
}

void CZIImage::closeImage()
{
  if( czi_reader != NULL ){
    czi_reader->Close();
    czi_reader = NULL;
  }
  if( tile_buf != NULL ){
    _TIFFfree( tile_buf );
    tile_buf = NULL;
  }
}



void CZIImage::loadImageInfo( int seq, int ang )
{
  currentX = seq;
  currentY = ang;

  auto mds = czi_reader->ReadMetadataSegment();
  auto md = mds->CreateMetaFromMetadataSegment();

  auto sbStatistics = czi_reader->GetStatistics();

  /// The number of channels for this image
  ///   IIPImage:  unsigned int channels;
  ///   CZIImage:  int channels_start;
  ///              int z_layers_start;
  ///              int z_layers_size;
  // See:  CZICmd/execute.cpp -- PrintStatistics().
  int channels_size = 0;  // Unchanged if TryGetInterval() == false.
  sbStatistics.dimBounds.TryGetInterval(libCZI::DimensionIndex::C, &channels_start, &channels_size);
  channels = (unsigned int) channels_size;

  sbStatistics.dimBounds.TryGetInterval(libCZI::DimensionIndex::Z, &z_layers_start, &z_layers_size);

  // Store the list of image dimensions available

  /// The image pixel dimensions
  ///   IIPImage:  std::vector <unsigned int> image_widths, image_heights;
  ///   CZIImage:  std::uint8_t image_minification;
  ///              std::vector <unsigned int> image_scales;
  /// The number of available resolutions in this image
  ///   IIPImage:  unsigned int numResolutions;

  // See:  CZICmd/execute.cpp -- PrintStatistics().
  /*
Bounding-Box:
 All:    X=0 Y=0 W=98172 H=66096
 Layer0: X=0 Y=0 W=98097 H=65963
  */
  unsigned int full_width = sbStatistics.boundingBox.w;
  unsigned int full_height = sbStatistics.boundingBox.h;
  ////unsigned int full_width = sbStatistics.boundingBoxLayer0Only.w;
  ////unsigned int full_height = sbStatistics.boundingBoxLayer0Only.h;
  image_widths.push_back( full_width );
  image_heights.push_back( full_height );
  image_minification = 2;  // Default to be updated.
  image_scales.push_back( 1 );
  numResolutions = 1;

  // CZI pyramid layers as virtual width x height images.
  // See:  CZICmd/execute.cpp -- PrintPyramidStatistics().
  /*
scene#0:
 number of subblocks with scale 1/1: 2550
 number of subblocks with scale 1/3: 919
 number of subblocks with scale 1/9: 116
 number of subblocks with scale 1/27: 17
 number of subblocks with scale 1/81: 3
 number of subblocks with scale 1/243: 1
  */
  auto pyrStatistics = czi_reader->GetPyramidStatistics();
  auto scene0PyrStatistics = pyrStatistics.scenePyramidStatistics[0];
  for (const auto& layer_stats : scene0PyrStatistics) {
    if (!layer_stats.layerInfo.IsNotIdentifiedAsPyramidLayer()) {
      if (layer_stats.layerInfo.IsLayer0() != true) {
        int scale = layer_stats.layerInfo.minificationFactor;
        for (int n = 1; n < layer_stats.layerInfo.pyramidLayerNo; ++n) {
          scale *= layer_stats.layerInfo.minificationFactor;
        }

        // Scaled down dimensions [rounded up with (.. -1)/scale +1]
        unsigned int w = (full_width -1)/scale +1;
        unsigned int h = (full_height -1)/scale +1;
        // Ignore downsamples smaller than 2K x 2K.
        // TODO(Leo) Why?  Ask Coleman.  Maybe don't care about viewing anything smaller.
        if (w < 2000 && h < 2000) {
          break;
        }

        image_widths.push_back( w );
        image_heights.push_back( h );
        image_minification = layer_stats.layerInfo.minificationFactor;
        image_scales.push_back( scale );
        ++numResolutions;
      }
    }
  }


  /// The base tile pixel dimensions
  ///   IIPImage:  unsigned int tile_width, tile_height;
  // CZI images are not stored as grids of tiles,
  // so output a default virtual tile dimensions, e.g.: 728x728 or 1024x1024
  // Virtual tiles will be created on the fly.
  tile_width = CZIIMAGE_DEFAULT_TILE_WIDTH;
  tile_height = CZIIMAGE_DEFAULT_TILE_HEIGHT;

  // See:  CZICmd/execute.cpp -- PrintAllSubBlocks().
  /*
Complete list of sub-blocks
---------------------------

#0: C0S0B0 M=0 logical=(0,0,98172,66096) phys.=(404,272) pixeltype=bgr24 comp.mode=jpgxr
  */
  libCZI::PixelType pixel_type = libCZI::PixelType::Invalid;
  czi_reader->EnumerateSubBlocks(
      [&](int index, const libCZI::SubBlockInfo& info)->bool
  {
    if (pixel_type == libCZI::PixelType::Invalid)
      pixel_type = info.pixelType;
    else if (pixel_type != info.pixelType)
      return false;

    return true;
  });

  /// The bits per channel for this image
  ///   unsigned int bpc;
  bpc = 8 * CziUtils::GetBytesPerPel(pixel_type);

  /// The colour space of the image
  ///   ColourSpaces colourspace;
  switch (pixel_type) {
    case libCZI::PixelType::Gray8:
    case libCZI::PixelType::Gray16:
    case libCZI::PixelType::Gray32:

    case libCZI::PixelType::Gray32Float:
    case libCZI::PixelType::Gray64Float:
    case libCZI::PixelType::Gray64ComplexFloat:
      colourspace = GREYSCALE;
      break;

    case libCZI::PixelType::Bgr24:
    case libCZI::PixelType::Bgr48:

    case libCZI::PixelType::Bgra32:

    case libCZI::PixelType::Bgr96Float:
    case libCZI::PixelType::Bgr192ComplexFloat:
      colourspace = sRGB;
      break;

    case libCZI::PixelType::Invalid:
    default:
      colourspace = NONE;
      break;
  }

  /// The sample format type (fixed or floating point)
  ///   SampleType sampleType;
  switch (pixel_type) {

    case libCZI::PixelType::Gray32Float:
    case libCZI::PixelType::Gray64Float:
    case libCZI::PixelType::Gray64ComplexFloat:

    case libCZI::PixelType::Bgr96Float:
    case libCZI::PixelType::Bgr192ComplexFloat:
      sampleType == FLOATINGPOINT;
      break;

    case libCZI::PixelType::Gray8:
    case libCZI::PixelType::Gray16:
    case libCZI::PixelType::Gray32:

    case libCZI::PixelType::Bgr24:
    case libCZI::PixelType::Bgr48:

    case libCZI::PixelType::Bgra32:

    case libCZI::PixelType::Invalid:
    default:
      sampleType = FIXEDPOINT;
      break;
  }

#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// The min and max sample value for each channel
  std::vector <float> min, max;  // Maybe --info-level 'DisplaySettings' ???
  /*
$ CZIcmd --command PrintInformation --info-level DisplaySettingsJSON --source ../18_WT1_NKPP002_01_R_LS.czi 
Display-Settings in CZIcmd-JSON-Format
--------------------------------------


Pretty-Print:
{
    "channels": [
        {
            "ch": 0,
            "black-point": 0.0,
            "white-point": 1.0
        }
    ]
}

Compact:
{"channels":[{"ch":0,"black-point":0.0,"white-point":1.0}]}


$ CZIcmd --command PrintInformation --info-level DisplaySettingsJSON --source ../GR57-13\ 2015_10_12__0078.czi 
Display-Settings in CZIcmd-JSON-Format
--------------------------------------


Pretty-Print:
{
    "channels": [
        {
            "ch": 0,
            "black-point": 0.0,
            "white-point": 0.900952398777008,
            "tinting": "#0000ff"
        },
        {
            "ch": 1,
            "black-point": 0.000053032286814413965,
            "white-point": 0.2518841624259949,
            "tinting": "#00ff5b"
        },
        {
            "ch": 2,
            "black-point": 0.00008445332059636712,
            "white-point": 0.2818259298801422,
            "tinting": "#ff0000"
        }
    ]
}

Compact:
{"channels":[{"ch":0,"black-point":0.0,"white-point":0.9009523987770081,"tinting":"#0000ff"},{"ch":1,"black-point":0.000053032286814413965,"white-point":0.2518841624259949,"tinting":"#00ff5b"},{"ch":2,"black-point":0.00008445332059636712,"white-point":0.2818259298801422,"tinting":"#ff0000"}]}
  */

#endif // TODO(Leo)  IIPImage metadata to be initialized.

#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// Quality layers
  //  unsigned int quality_layers;  // Not used for CZI.

  /// Image histogram
  //  std::vector<unsigned int> histogram;  // Not needed here.

#endif // TODO(Leo)  IIPImage metadata to be initialized.

#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// STL map to hold string metadata
  std::map <const std::string, std::string> metadata;
  /*
$ CZIcmd --command PrintInformation --info-level GeneralInfo --source ../GR57-13\ 2015_10_12__0078.czi 
General Information
-------------------

Name=
Title=
UserName=cgib
Description=
Comment=
Keywords=
Rating=0
CreationDate=2015-10-12T14:05:54.4916631-07:00
  */
#endif // TODO(Leo)
}


// See:  CZICmd/SaveBitmap.cpp -- CSaveData::SaveBgr24()
void CZIImage::tweakLine(libCZI::PixelType pixel_type, std::uint32_t width, void* ptrData) {
  if (pixel_type == libCZI::PixelType::Bgr24) {
	char* p = (char*)ptrData;
	for (std::uint32_t x = 0; x < width; ++x) {
	  std::swap(*p,*(p+2));
	  p+=3;
	}
  }
}

RawTile CZIImage::getTile( int seq, int ang, unsigned int res, int layers, unsigned int tile )
{
  logfile << "CZIImage:: getTile begin" << endl;

  // Check the resolution exists
  if( res > (numResolutions - 1) ){
    ostringstream error;
    error << "CZIImage :: Asked for non-existent resolution: " << res;
    throw file_error( error.str() );
  }

  // If we are currently working on a different sequence number,
  // then close and reload the image.
  if( (currentX != seq) || (currentY != ang) ){
    closeImage();
  }

  // Allocate and open the CZI if it's not already open.
  if( !czi_reader ){
    string filename = getFileName( seq, ang );
    updateTimestamp( filename );
    std::wstring wide_filename = std::wstring(filename.begin(), filename.end());
    auto stream = libCZI::CreateStreamFromFile(wide_filename.c_str());
    czi_reader = libCZI::CreateCZIReader();
    czi_reader->Open(stream);
  }

  // Reload our image information in case the tile size etc is different
  if( (currentX != seq) || (currentY != ang) ){
    loadImageInfo( seq, ang );
  }


  // The first CZI layer is the highest resolution,
  // so we need to invert the resolution index:
  //   0 [lowest resolution] <= res <= (numResolutions - 1) [highest resolution]
  int czi_layer = (numResolutions - 1) - res;

  unsigned int image_width = image_widths[czi_layer];
  unsigned int image_height = image_heights[czi_layer];

  // CZI images are not stored as grids of tiles,
  // so create virtual tile grid on the fly.

#if 0  // TODO(Leo)
  /*
                CZIcmd --source "${UPLOAD}" \
                    --command ChannelComposite \
                    --rect "rel($X,$Y,$W,$H)" \
                    ${PLANE} ${DISPLAY_CHANNEL} \
                    --background ${BACKGROUND} \
                    --output "${OUTPUT}"

  /// For brightfield images, using SingleChannelPyramidTileAccessor works well.
CZIcmd --source ../../8_PPIB_NKPP002_01_R_LS.czi \
    --command SingleChannelPyramidTileAccessor \
    --rect 'rel(40960,40960,2048,2048)' \
    --background 0.9 \
    --plane-coordinate C0 \
    --pyramidinfo 3,0 \
    --output 8_PPIB_NKPP002_01_R_LS.czi--command-SingleChannelPyramidTileAccessor--rect-rel40960,40960,2048,2048--background0.9--plane-coordinate-C0--pyramidinfo-3,0

CZIcmd --source ../../8_PPIB_NKPP002_01_R_LS.czi \
    --command SingleChannelPyramidTileAccessor \
    --rect 'rel(40960,40960,6144,6144)' \
    --background 0.9 \
    --plane-coordinate C0 \
    --pyramidinfo 3,1 \
    --output 8_PPIB_NKPP002_01_R_LS.czi--command-SingleChannelPyramidTileAccessor--rect-rel40960,40960,6144,6144--background0.9--plane-coordinate-C0--pyramidinfo-3,1

CZIcmd --source ../../8_PPIB_NKPP002_01_R_LS.czi \
    --command SingleChannelPyramidTileAccessor \
    --rect 'rel(40960,40960,18432,18432)' \
    --background 0.9 \
    --plane-coordinate C0 \
    --pyramidinfo 3,2 \
    --output 8_PPIB_NKPP002_01_R_LS.czi--command-SingleChannelPyramidTileAccessor--rect-rel40960,40960,18432,18432--background0.9--plane-coordinate-C0--pyramidinfo-3,2

CZIcmd --source ../../8_PPIB_NKPP002_01_R_LS.czi \
    --command SingleChannelPyramidTileAccessor \
    --rect 'rel(0,0,102546,73629)' \
    --background 0.9 \
    --plane-coordinate C0 \
    --pyramidinfo 3,4 \
    --output 8_PPIB_NKPP002_01_R_LS.czi--command-SingleChannelPyramidTileAccessor--rect-rel0,0,102546,73629--background0.9--plane-coordinate-C0--pyramidinfo-3,4


    /// For fluorescence images, ScalingChannelComposite works well.
CZIcmd -s../../GR57-13\ 2015_10_12__0078.czi \
    -cScalingChannelComposite \
    -r'rel(20480,20480,32768,32768)' \
    -b0.0 \
    -pC0 \
    -d'{"channels":[{"ch":0,"black-point":0.0,"white-point":0.9009523987770081,"tinting":"#ffffff"}]}'  \
    -z 0.0625 \
    --output GR57-13\ 2015_10_12__0078.czi-cScalingChannelComposite-rrel20480,20480,32768,32768-b0.0-pC0-dchannels..ch0..ffffff-z0.0625
  */
#endif // TODO(Leo)

  // [Use (.. -1)/.. +1 to round up the division.]
  unsigned int num_cols = (image_width -1)/tile_width +1;
  unsigned int num_rows = (image_height -1)/tile_height +1;

  // Check that a valid tile number was given
  if( tile >= num_cols * num_rows ) {
    ostringstream tile_no;
    tile_no << "Asked for non-existent tile: " << tile;
    throw file_error( tile_no.str() );
  }

  // Get grid (col, row) and tile upper-left (x, y) coordinates.
  unsigned int grid_col = tile % num_cols;
  unsigned int grid_row = tile / num_cols;
  unsigned int tile_x = grid_col * tile_width;
  unsigned int tile_y = grid_row * tile_height;

  // Get pixel dimensions of this tile, which can be smaller than the
  // default tile dimensions for edge tiles.
  unsigned int tile_w = image_width - tile_x;  // Distance to edge.
  if (tile_width < tile_w) tile_w = tile_width;
  unsigned int tile_h = image_height - tile_y;  // Distance to edge.
  if (tile_height < tile_h) tile_h = tile_height;

  // Number of pixels for this tile.
  unsigned int num_pixels = tile_w * tile_h; 


  /// For brightfield images, using SingleChannelPyramidTileAccessor works well.
  // See:  CZICmd/execute.cpp -- class CExecuteSingleChannelPyramidTileAccessor -- execute().
  auto subBlockStatistics = czi_reader->GetStatistics();
  auto accessor = czi_reader->CreateSingleChannelPyramidLayerTileAccessor();

  unsigned int scale = image_scales[czi_layer];
  int logical_x = tile_x * scale + subBlockStatistics.boundingBox.x;
  int logical_y = tile_y * scale + subBlockStatistics.boundingBox.y;
  int logical_w = tile_w * scale;
  int logical_h = tile_h * scale;
  libCZI::IntRect roi{logical_x, logical_y, logical_w, logical_h};

  libCZI::CDimCoordinate coordinate;
  if (channels > 0)
    coordinate.Set(libCZI::DimensionIndex::C, channels_start + 0);  // For now, just first Channel.
  if (z_layers_size > 0)
    coordinate.Set(libCZI::DimensionIndex::Z, z_layers_start + 0);  // For now, just first Z-layer.

  libCZI::ISingleChannelPyramidLayerTileAccessor::Options scptaOptions;
  scptaOptions.Clear();
  libCZI::RgbFloatColor bright_bkgd{ 0.9, 0.9, 0.9 };  // Light for brightfield images.
  libCZI::RgbFloatColor fluor_bkgd{ 0.0, 0.0, 0.0 };  // Black for fluorescence channels.
  scptaOptions.backGroundColor = (channels == 1 ? bright_bkgd : fluor_bkgd);
  //  scptaOptions.sceneFilter = options.GetSceneIndexSet(); // Unused, leave as default.

  libCZI::ISingleChannelPyramidLayerTileAccessor::PyramidLayerInfo pyrLyrInfo;
  pyrLyrInfo.minificationFactor = image_minification;
  pyrLyrInfo.pyramidLayerNo = czi_layer;

  // std::shared_ptr<libCZI::IBitmapData>
  auto bitmap = accessor->Get(roi, &coordinate, pyrLyrInfo, &scptaOptions);

  // See:  CZICmd/SaveBitmap.cpp
  //    -- CSaveData::Save(), CSaveData::SaveBgr24(), CSaveData::SaveGray16(),
  //    -- CSaveData::SavePngTweakLineBeforeWritng(), CSaveData::SavePng()
  libCZI::ScopedBitmapLockerP lckScoped{bitmap.get()};

  // Allocate memory for our tile.
  if( !tile_buf ){
    if( ( tile_buf = _TIFFmalloc( lckScoped.size ) ) == NULL ){
      throw file_error( "tiff malloc tile failed" );
    }
  }

  std::unique_ptr<void,decltype(&free)> lineToTweak(malloc(lckScoped.stride),&free);
  for (std::uint32_t h = 0; h < bitmap->GetHeight(); ++h) {
    void *roi_ptr = (((char *) lckScoped.ptrDataRoi) + h * lckScoped.stride);
    memcpy(lineToTweak.get(), roi_ptr, lckScoped.stride);
	tweakLine(bitmap->GetPixelType(), bitmap->GetWidth(), lineToTweak.get());

    void *tile_buf_ptr = (((char *) tile_buf) + h * lckScoped.stride);
    memcpy(tile_buf_ptr, roi_ptr, lckScoped.stride);
  }

  RawTile rawtile( tile, res, seq, ang, tile_w, tile_h, channels, bpc );
  rawtile.data = tile_buf;
  rawtile.dataLength = lckScoped.size;
  rawtile.filename = getImagePath();
  rawtile.timestamp = timestamp;
  rawtile.memoryManaged = 0;
  //rawtile.padded = true;  // TODO(Leo) Huh? Why padded for QPTIFF, not for OpenSlide?
  rawtile.sampleType = sampleType;







  /// For fluorescence images, ScalingChannelComposite works well.
  
  return rawtile;
}
