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
#include "CziUtils.h"


#if 1  // For debugging to logfile (/tmp/iipsrv.log).
#include <fstream>  // operator<<(), __FILE__,...
extern int loglevel;
extern std::ofstream logfile;
#endif // For debugging to logfile (/tmp/iipsrv.log).

using namespace std;  // string, endl


void CZIImage::openImage() {
  if (loglevel >= 3)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  BEGIN" << endl;

  // Insist that the czi_reader and tile_buf be NULL.
  if ( czi_reader || tile_buf ) {
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << ", czi_reader or tile_buf is not NULL" << endl;
    throw file_error( "czi_reader or tile_buf is not NULL" );
  }

  string filename = getFileName( currentX, currentY );

  // Update our timestamp.
  updateTimestamp( filename );

  // Try to allocate and open a CZI-reader.
  // See:  libCZI/Src/CZICmd/execute.cpp -- class CExecuteBase -- CreateAndOpenCziReader().
  std::wstring wide_filename = std::wstring(filename.begin(), filename.end());
  auto stream = libCZI::CreateStreamFromFile(wide_filename.c_str());
  czi_reader = libCZI::CreateCZIReader();
  czi_reader->Open(stream);

  // Load our metadata if not already loaded.
  if ( bpc == 0 || image_scales.size() < 1 ) {
    loadImageInfo( currentX, currentY );
  }

  // Insist on a tiled image.
  if ( (tile_width == 0) && (tile_height == 0) ) {
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << ", CZI image is not tiled" << endl;
    throw file_error( "CZI image is not tiled" );
  }

  isSet = true;

  if (loglevel >= 5)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  END" << endl;
}

void CZIImage::closeImage() {
  if (loglevel >= 3)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  BEGIN" << endl;

  if ( czi_reader != NULL ) {
    czi_reader->Close();
    czi_reader = NULL;
  }
  if ( tile_buf != NULL ) {
    _TIFFfree( tile_buf );
    tile_buf = NULL;
    tile_buf_size = 0;
  }

  if (loglevel >= 5)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  END" << endl;
}

void CZIImage::tile_buf_malloc(tmsize_t size) {
  if (loglevel >= 11)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "( size " << size << " )  BEGIN" << endl;

  // If increasing the size of the tile buffer, delete the old one.
  if (tile_buf && tile_buf_size < size) {
    _TIFFfree( tile_buf );
    tile_buf = NULL;
    tile_buf_size = 0;
  }

  // Allocate memory for our tile.
  if ( !tile_buf ) {
    if ( ( tile_buf = _TIFFmalloc( size ) ) == NULL ) {
      throw file_error( "tiff malloc tile failed" );
    }
    tile_buf_size = size;
  }

  if (loglevel >= 11)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  END" << endl;
}


void CZIImage::loadImageInfo( int seq, int ang ) {
  if (loglevel >= 3)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  BEGIN" << endl;

  currentX = seq;
  currentY = ang;


  auto subblock_statistics = czi_reader->GetStatistics();

  // Store the list of image dimensions available

  /// The image pixel dimensions
  ///   IIPImage:  std::vector <unsigned int> image_widths, image_heights;
  ///   CZIImage:  std::uint8_t image_minification;
  ///              std::vector <unsigned int> image_scales;
  /// The number of available resolutions in this image
  ///   IIPImage:  unsigned int numResolutions;
  //
  // See:  libCZI/Src/CZICmd/execute.cpp -- PrintStatistics().
  /*
Bounding-Box:
 All:    X=0 Y=0 W=98172 H=66096
 Layer0: X=0 Y=0 W=98097 H=65963
  */
  unsigned int full_width = subblock_statistics.boundingBox.w;
  unsigned int full_height = subblock_statistics.boundingBox.h;
  ////unsigned int full_width = subblock_statistics.boundingBoxLayer0Only.w;
  ////unsigned int full_height = subblock_statistics.boundingBoxLayer0Only.h;

  image_widths.push_back( full_width );
  image_heights.push_back( full_height );
  image_minification = 2;  // Default to be updated.
  image_scales.push_back( 1 );
  numResolutions = 1;

  if (loglevel >= 3)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << ", full_width " << full_width << ", full_height " << full_height
            << endl;

  // CZI pyramid layers as virtual width x height images.
  // (Actually stored as CZI subblocks - arbitrarily aligned scaled tiles.)
  // See:  libCZI/Src/CZICmd/execute.cpp -- PrintPyramidStatistics().
  /*
scene#0:
 number of subblocks with scale 1/1: 2550
 number of subblocks with scale 1/3: 919
 number of subblocks with scale 1/9: 116
 number of subblocks with scale 1/27: 17
 number of subblocks with scale 1/81: 3
 number of subblocks with scale 1/243: 1
  */
  /*
scene#1:
 number of subblocks with scale 1/1: 114
 number of subblocks with scale 1/2: 54
 number of subblocks with scale 1/4: 16
 number of subblocks with scale 1/8: 4
 number of subblocks with scale 1/16: 1

scene#2:
 number of subblocks with scale 1/1: 1306
 number of subblocks with scale 1/2: 539
 number of subblocks with scale 1/4: 149
 number of subblocks with scale 1/8: 44
 number of subblocks with scale 1/16: 14
 number of subblocks with scale 1/32: 4
 number of subblocks with scale 1/64: 1

scene#3:
 number of subblocks with scale 1/1: 2084
 number of subblocks with scale 1/2: 876
 number of subblocks with scale 1/4: 237
 number of subblocks with scale 1/8: 69
 number of subblocks with scale 1/16: 20
 number of subblocks with scale 1/32: 7
 number of subblocks with scale 1/64: 3
 number of subblocks with scale 1/128: 1
  */
  auto pyramid_statistics = czi_reader->GetPyramidStatistics();

  // Find the largest common scaling denominator used by all subblocks.
  // See:  libCZI/Src/CZICmd/execute.cpp -- PrintPyramidStatistics().
  int common_scale_denominator = -1;
  for (const auto& i : pyramid_statistics.scenePyramidStatistics) {
    // scene-index==int::max means "scene-index not valid"
    if (i.first != (numeric_limits<int>::max)()) {
      if (loglevel >= 3)
        logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
                << "scene#" << i.first << ":" << endl;
    }

    int scale_denominator;
    for (const auto& j : i.second) {
      if (!j.layerInfo.IsNotIdentifiedAsPyramidLayer()) {
        if (j.layerInfo.IsLayer0() == true) {
          scale_denominator = 1;
        }
        else {
          scale_denominator = j.layerInfo.minificationFactor;
          for (int n = 0; n < j.layerInfo.pyramidLayerNo - 1; ++n)
            {
              scale_denominator *= j.layerInfo.minificationFactor;
            }
        }

        if (loglevel >= 3)
          logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
                  << " number of subblocks with scale 1/" << scale_denominator << ": " << j.count << endl;
      }
      else {
        if (loglevel >= 3)
          logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
                  << " number of subblocks not representable as pyramid-layers: " << j.count << endl;;
      }
    }
    if (common_scale_denominator < 0 || common_scale_denominator > scale_denominator) {
      common_scale_denominator = scale_denominator;
    }
  }
  if (loglevel >= 3)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << ", common_scale_denominator " << common_scale_denominator
            << endl;

  auto scene0_pyramid_statistics = pyramid_statistics.scenePyramidStatistics[0];
  for (const auto& layer_statistics : scene0_pyramid_statistics) {
    if (layer_statistics.layerInfo.IsNotIdentifiedAsPyramidLayer() != true) {
      if (layer_statistics.layerInfo.IsLayer0() != true) {
        int scale = layer_statistics.layerInfo.minificationFactor;
        for (int n = 1; n < layer_statistics.layerInfo.pyramidLayerNo; ++n) {
          scale *= layer_statistics.layerInfo.minificationFactor;
        }
        if (scale > common_scale_denominator) {
          if (loglevel >= 3)
            logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
                    << ", skipping scale " << scale
                    << endl;

          break;
        }

        /*TEMP(Leo)*//* if (layer_statistics.layerInfo.minificationFactor != 2) { break; } */

        // Scaled down dimensions [rounded up with (.. -1)/scale +1]
        unsigned int w = (full_width -1)/scale +1;
        unsigned int h = (full_height -1)/scale +1;

#if 1  // It's a QPTIFF thing, smaller than 2K x 2K are strips, not tiles.
        // Ignore downsamples smaller than 2K x 2K.
        // TODO(Leo) Why?  Ask Coleman.  Maybe don't care about viewing anything smaller.
        if (w < 2000 && h < 2000) {
          break;
        }
#endif // It's a QPTIFF thing, smaller than 2K x 2K are strips, not tiles.

        // High resolution to low resolution.
        image_widths.push_back( w );
        image_heights.push_back( h );
        image_minification = layer_statistics.layerInfo.minificationFactor;
        resolution_scale_factor = image_minification;
        image_scales.push_back( scale );
        ++numResolutions;

        if (loglevel >= 3)
          logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
                  << ", w " << w << ", h " << h
                  << ", image_minification " << (int)image_minification
                  << ", numResolutions " << numResolutions
                  << endl;
      }
    }
  }

  if (loglevel >= 3)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << ", image_minification " << (int)image_minification
            << ", numResolutions " << numResolutions
            << endl;


  /// The base tile pixel dimensions
  ///   IIPImage:  unsigned int tile_width, tile_height;
  //
  // CZI images are not stored as grids of tiles,
  // so output a default virtual tile dimensions, e.g.: 728x728 or 1024x1024
  // Virtual tiles will be created on the fly.
  // (Actually stored as CZI subblocks - arbitrarily aligned scaled tiles.)
  tile_width = CZIIMAGE_DEFAULT_TILE_WIDTH;
  tile_height = CZIIMAGE_DEFAULT_TILE_HEIGHT;


  // See:  libCZI/Src/CZICmd/execute.cpp -- PrintAllSubBlocks().
  /*
Complete list of sub-blocks
---------------------------

#0: C0S0B0 M=0 logical=(0,0,98172,66096) phys.=(404,272) pixeltype=bgr24 comp.mode=jpgxr
  */
  libCZI::PixelType pixel_type = libCZI::PixelType::Invalid;
  czi_reader->EnumerateSubBlocks([&](int index, const libCZI::SubBlockInfo& info)->bool
                                 {
                                   if (pixel_type == libCZI::PixelType::Invalid)
                                     pixel_type = info.pixelType;
                                   else if (pixel_type != info.pixelType)
                                     return false;

                                   return true;
                                 });

  /// The colour space of the image
  ///   ColourSpaces colourspace;
  //
  // See:  libCZI/Src/libCZI/libCZI_Pixels.h -- enum class PixelType.
  switch (pixel_type) {
    case libCZI::PixelType::Gray8:  case libCZI::PixelType::Gray16:  case libCZI::PixelType::Gray32:
    case libCZI::PixelType::Gray32Float:  case libCZI::PixelType::Gray64Float:  case libCZI::PixelType::Gray64ComplexFloat:
      colourspace = GREYSCALE;
      break;

    case libCZI::PixelType::Bgr24:  case libCZI::PixelType::Bgr48:  case libCZI::PixelType::Bgra32:
    case libCZI::PixelType::Bgr96Float:  case libCZI::PixelType::Bgr192ComplexFloat:
      colourspace = sRGB;
      break;

    case libCZI::PixelType::Invalid:  default:
      colourspace = NONE;
      break;
  }

  /// The sample format type (fixed or floating point)
  ///   SampleType sampleType;
  switch (pixel_type) {
    case libCZI::PixelType::Gray32Float:  case libCZI::PixelType::Gray64Float:  case libCZI::PixelType::Gray64ComplexFloat:
    case libCZI::PixelType::Bgr96Float:  case libCZI::PixelType::Bgr192ComplexFloat:
      sampleType == FLOATINGPOINT;
      break;

    case libCZI::PixelType::Gray8:  case libCZI::PixelType::Gray16:  case libCZI::PixelType::Gray32:
    case libCZI::PixelType::Bgr24:  case libCZI::PixelType::Bgr48:  case libCZI::PixelType::Bgra32:
    case libCZI::PixelType::Invalid:  default:
      sampleType = FIXEDPOINT;
      break;
  }

#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// The min and max sample value for each channel
  std::vector <float> min, max;  // Maybe --info-level 'DisplaySettings' ???
  /*
$ CZIcmd --command PrintInformation --info-level DisplaySettingsJSON --source ../GR57-13\ 2015_10_12__0078.czi 
Display-Settings in CZIcmd-JSON-Format
...
{"channels":[{"ch":0,"black-point":0.0,"white-point":0.9009523987770081,"tinting":"#0000ff"},{"ch":1,"black-point":0.000053032286814413965,"white-point":0.2518841624259949,"tinting":"#00ff5b"},{"ch":2,"black-point":0.00008445332059636712,"white-point":0.2818259298801422,"tinting":"#ff0000"}]}
  */
#endif // TODO(Leo)  IIPImage metadata to be initialized.

#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// Quality layers
  unsigned int quality_layers;  // Not used for CZI.

  /// Image histogram
  std::vector<unsigned int> histogram;  // Not needed here.
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


  /// The CZI dimensions for this image
  ///   CZIImage:  int channels_start;
  ///              int channels_size;
  ///              int z_layers_start;
  ///              int z_layers_size;
  //
  // See:  libCZI/Src/CZICmd/execute.cpp -- PrintStatistics().
  subblock_statistics.dimBounds.TryGetInterval(libCZI::DimensionIndex::C,
                                               &channels_start, &channels_size);
  subblock_statistics.dimBounds.TryGetInterval(libCZI::DimensionIndex::Z,
                                               &z_layers_start, &z_layers_size);
  if (loglevel >= 3)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << ", channels_start " << channels_start << ", channels_size " << channels_size
            << ", z_layers_start " << z_layers_start << ", z_layers_size " << z_layers_size
            << endl;


  /// The number of IIPImage channels for this image
  ///   IIPImage:  unsigned int channels;
  channels = (unsigned int) (channels_size < 1 ? 1 : channels_size);  // Handle undefined C dimension.

//#define USE_ONLY_ACTIVE_ENABLED_CHANNELS
#undef USE_ONLY_ACTIVE_ENABLED_CHANNELS  // #undef to use all channels
#ifdef USE_ONLY_ACTIVE_ENABLED_CHANNELS
  if (channels > 1) {
    std::shared_ptr<libCZI::IDisplaySettings> full_display_settings = (czi_reader
                                                                       ->ReadMetadataSegment()
                                                                       ->CreateMetaFromMetadataSegment()
                                                                       ->GetDocumentInfo()
                                                                       ->GetDisplaySettings());
    std::vector<int> active_channels = libCZI::CDisplaySettingsHelper::GetActiveChannels(full_display_settings.get());
    channels = active_channels.size();
  }
#endif

  /// The bits per channel for this image
  ///   unsigned int bpc;
  bpc = 8 * CziUtils::GetBytesPerPel(pixel_type);


  // Clean-up:  libCZI and IIPImage terminology about channels and bits per sample disagree.
  //            Also, "JPEGCompressor: JPEG can only handle 8 bit images"
  if (loglevel >= 2)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << ", CZI Reader: channels " << channels << ", bpc " << bpc << endl;

  if (pixel_type == libCZI::PixelType::Bgr24) {
    // Single CZI Brg24 channel to be rendered as collated IIPSRV 8-bit R/G/B channels.
    channels *= (bpc / 8);
    bpc = 8;
  }

  if (loglevel >= 2)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << ", IIP Server: channels " << channels << ", bpc " << bpc << endl;

  if (loglevel >= 5)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  END" << endl;
}


// CZI images are not stored as grids of tiles,
// so create virtual tile grid on the fly.
//
// Samples of creating tiles of regions using CZIcmd:
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
    --output "${OUTPUT}"

CZIcmd --source ../../8_PPIB_NKPP002_01_R_LS.czi \
    --command SingleChannelPyramidTileAccessor \
    --rect 'rel(40960,40960,6144,6144)' \
    --background 0.9 \
    --plane-coordinate C0 \
    --pyramidinfo 3,1 \
    --output "${OUTPUT}"

CZIcmd --source ../../8_PPIB_NKPP002_01_R_LS.czi \
    --command SingleChannelPyramidTileAccessor \
    --rect 'rel(40960,40960,18432,18432)' \
    --background 0.9 \
    --plane-coordinate C0 \
    --pyramidinfo 3,2 \
    --output "${OUTPUT}"

CZIcmd --source ../../8_PPIB_NKPP002_01_R_LS.czi \
    --command SingleChannelPyramidTileAccessor \
    --rect 'rel(0,0,102546,73629)' \
    --background 0.9 \
    --plane-coordinate C0 \
    --pyramidinfo 3,4 \
    --output "${OUTPUT}"


/// For fluorescence images, ScalingChannelComposite works well.
CZIcmd -s../../GR57-13\ 2015_10_12__0078.czi \
    --command ScalingChannelComposite \
    --rect 'rel(20480,20480,32768,32768)' \
    --background 0.0 \
    --plane-coordinate C0 \
    --display-settings '{"channels":[{"ch":0,"black-point":0.0,"white-point":0.9009523987770081,"tinting":"#ffffff"}]}'  \
    --zoom 0.0625 \
    --output "${OUTPUT}"
*/

RawTile CZIImage::getTile( int seq, int ang, unsigned int res, int unused_layers, unsigned int tile )
{
  if (loglevel >= 5)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__
            << "( seq " << seq << ", ang " << ang << ", res " << res
            << ", layers " << unused_layers << ", tile " << tile << " )  BEGIN" << endl;

  // Check the resolution exists
  if ( res > (numResolutions - 1) ) {
    ostringstream res_error;
    res_error << "Asked for non-existent resolution: " << res;
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "(), " << res_error.str() << endl;
    throw file_error( res_error.str() );
  }

  // If we are currently working on a different sequence number,
  // then close and reload the image.
  if ( (currentX != seq) || (currentY != ang) ) {
    closeImage();
  }

  // Allocate and open the CZI if it's not already open.
  if ( !czi_reader ) {
    currentX = seq;
    currentY = ang;
    bpc = 0;
    openImage();
  }

  // Reload our image information in case the tile size, etc. is different.
  if ( (currentX != seq) || (currentY != ang) || image_scales.size() < 1 ) {
    loadImageInfo( seq, ang );
  }


  // The first CZI layer is the highest resolution,
  // so we need to invert the IIPSRV resolution index:
  //   0 [lowest resolution] <= res <= (numResolutions - 1) [highest resolution]
  int czi_pyr_layer = (numResolutions - 1) - res;

  // Image and grid dimensions for selected resolution.
  unsigned int image_width = image_widths[czi_pyr_layer];
  unsigned int image_height = image_heights[czi_pyr_layer];
  // [Use (.. -1)/.. +1 to round up the division.]
  unsigned int num_cols = (image_width -1)/tile_width +1;
  unsigned int num_rows = (image_height -1)/tile_height +1;

  // Check that a valid tile number was given
  if ( tile >= num_cols * num_rows ) {
    ostringstream tile_error;
    tile_error << "Asked for non-existent tile: " << tile;
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "(), " << tile_error.str() << endl;
    throw file_error( tile_error.str() );
  }

  // Get grid (col, row) and tile upper-left (x, y) coordinates.
  unsigned int grid_col = tile % num_cols;
  unsigned int grid_row = tile / num_cols;
  unsigned int tile_x = grid_col * tile_width;
  unsigned int tile_y = grid_row * tile_height;

  // Get pixel dimensions of this tile, which can be smaller than the
  // default tile dimensions for edge tiles.
  unsigned int tile_w = image_width - tile_x;  // Distance to right edge.
  if (tile_width < tile_w) { tile_w = tile_width; }  // Default width, if smaller.
  unsigned int tile_h = image_height - tile_y;  // Distance to bottom edge.
  if (tile_height < tile_h) { tile_h = tile_height; }  // Default height, if smaller.

  // Number of pixels for this tile.
  unsigned int num_pixels = tile_w * tile_h; 


  unsigned int scale = image_scales[czi_pyr_layer];
  auto subblock_statistics = czi_reader->GetStatistics();

  // Pixel coordinates in original (highest resolution) image.
  int logical_x = tile_x * scale + subblock_statistics.boundingBox.x;
  int logical_y = tile_y * scale + subblock_statistics.boundingBox.y;
  int logical_w = tile_w * scale;
  int logical_h = tile_h * scale;
  libCZI::IntRect roi{logical_x, logical_y, logical_w, logical_h};

  if (loglevel >= 2)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << ", tile x,y,w,h (" << tile_x << "," << tile_y << "," << tile_w << "," << tile_h << ")"
            << ", scale " << scale << ", logical roi " << roi
            << endl;

  int z_layer = 0;  // For now, just first Z-layer.

  if (channels_size == 1 /* Undefined C dimension: */ || channels_size < 0) {
    // Single CZI channel (usu. brightfield, Bgr24).
    return getSingleChannelPyramidLayerTile(seq,  ang, res, tile,
                                            z_layer, czi_pyr_layer,
                                            tile_w, tile_h, roi);
  }
  else if (channels_size > 1 /* Undefined C dimension: || channels_size < 0 */) {
    // Multiple CZI channels (usu. fluorescence, Gray16).
    return getAllChannelsPyramidLayerTile(seq, ang, res, tile,
                                          z_layer, czi_pyr_layer,
                                          tile_w, tile_h, roi);
  }
  else {
    // Throw an error -- should never get here.
    ostringstream channels_error;
    channels_error << "Invalid channels size: " << channels_size;
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "(), " << channels_error.str() << endl;
    throw file_error( channels_error.str() );
  }

  if (loglevel >= 5)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  END" << endl;
  return RawTile();
}


// See:  libCZI/Src/CZICmd/SaveBitmap.cpp -- CSaveData::SaveBgr24()
void CZIImage::tweakLine(const libCZI::PixelType pixel_type, const std::uint32_t width, void* const ptrData) {
  if (loglevel >= 11)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()" << endl;
  if (pixel_type == libCZI::PixelType::Bgr24) {
    char* p = (char*)ptrData;
    for (std::uint32_t x = 0; x < width; ++x) {
      std::swap(*p,*(p+2));
      p+=3;
    }
  }
}

// Single CZI channel (usu. brightfield, Bgr24).
RawTile CZIImage::getSingleChannelPyramidLayerTile(
    const int seq, const int ang, const unsigned int res, const unsigned int tile,
    const int z_layer, const int czi_pyr_layer,
    const unsigned int tile_w, const unsigned int tile_h, const libCZI::IntRect roi)
{
  if (loglevel >= 2)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__
            << "( seq " << seq << ", ang " << ang << ", res " << res << ", tile " << tile
            << ", z_layer " << z_layer << ", czi_pyr_layer " << czi_pyr_layer
            << ", tile_w " << tile_w << ", tile_h " << tile_h << ", roi " << roi
            << " )  BEGIN" << endl;

  /// For brightfield images, using Command::SingleChannelPyramidTileAccessor works well.
  // See:  libCZI/Src/CZICmd/execute.cpp -- class CExecuteSingleChannelPyramidTileAccessor -- execute().

  libCZI::CDimCoordinate coordinate;
  if (channels_size > 0) {
    coordinate.Set(libCZI::DimensionIndex::C, channels_start);
  }
  if (z_layers_size > 0) {
    coordinate.Set(libCZI::DimensionIndex::Z, z_layers_start + z_layer);
  }

  libCZI::ISingleChannelPyramidLayerTileAccessor::PyramidLayerInfo pyr_layer_info;
  pyr_layer_info.minificationFactor = image_minification;
  pyr_layer_info.pyramidLayerNo = czi_pyr_layer;

  libCZI::ISingleChannelPyramidLayerTileAccessor::Options scpta_options;
  scpta_options.Clear();
  //  scpta_options.sceneFilter = options.GetSceneIndexSet(); // Unused, leave as default from Clear().
  libCZI::RgbFloatColor bright_bkgd{ 1, 1, 1 /*0.9, 0.9, 0.9*/ };  // Light for brightfield images.
  /*TEMP(for background debugging)  libCZI::RgbFloatColor bright_bkgd{ 0.7,0.1,0.7 };*/  // Light for brightfield images.
  libCZI::RgbFloatColor fluor_bkgd{ 0.0, 0.0, 0.0 };  // Black for fluorescence channels.
  scpta_options.backGroundColor = (channels_size == 1 ? bright_bkgd : fluor_bkgd);
  /*TEMP(for subblock debugging)  scpta_options.drawTileBorder = true;*/

  auto accessor = czi_reader->CreateSingleChannelPyramidLayerTileAccessor();
  // std::shared_ptr<libCZI::IBitmapData>
  auto bitmap = accessor->Get(roi, &coordinate, pyr_layer_info, &scpta_options);

  if (loglevel >= 2)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << "  bitmap: width " << bitmap->GetWidth() << ", height " << bitmap->GetHeight()
            << ", pixel type " << (int) bitmap->GetPixelType() << ", size " << bitmap->GetSize()
            << endl;

  // See:  libCZI/Src/CZICmd/SaveBitmap.cpp
  //    -- CSaveData::Save(), CSaveData::SaveBgr24(), CSaveData::SaveGray16(),
  //    -- CSaveData::SavePngTweakLineBeforeWritng(), CSaveData::SavePng()
  libCZI::ScopedBitmapLockerP locked_bitmap{bitmap.get()};
  std::uint32_t tile_buf_stride =
    (locked_bitmap.stride / bitmap->GetWidth()) * bitmap->GetWidth();

  tile_buf_malloc( locked_bitmap.size );

  if (loglevel >= 2)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << "  locked_bitmap: stride " << locked_bitmap.stride << ", size " << locked_bitmap.size
            << "; tile_buf: stride " << tile_buf_stride << ", size " << tile_buf_size
            << endl;

  std::unique_ptr<void, decltype(&free)> lineToTweak(malloc(locked_bitmap.stride), &free);
  for (std::uint32_t h = 0; h < bitmap->GetHeight(); ++h) {
    void *roi_ptr = (((char *) locked_bitmap.ptrDataRoi) + h * locked_bitmap.stride);
    memcpy(lineToTweak.get(), roi_ptr, locked_bitmap.stride);
    tweakLine(bitmap->GetPixelType(), bitmap->GetWidth(), lineToTweak.get());

    void *tile_buf_ptr = (((char *) tile_buf) + h * tile_buf_stride);
    memcpy(tile_buf_ptr, lineToTweak.get(), tile_buf_stride);
  }

  if (loglevel >= 2)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << "  RawTile(): tile " << tile << ", res " << res
            << ", seq " << seq << ", ang " << ang
            << ", tile_w " << tile_w << ", tile_h " << tile_h
            << ", channels " << channels << ", bpc " << bpc
            << endl;

  RawTile rawtile( tile, res, seq, ang, tile_w, tile_h, channels, bpc );
  rawtile.data = tile_buf;
  rawtile.dataLength = locked_bitmap.size;
  rawtile.filename = getImagePath();
  rawtile.timestamp = timestamp;
  rawtile.memoryManaged = 0;
  //rawtile.padded = true;  // TODO(Leo) Huh? Why padded for QPTIFF, not for OpenSlide?
  rawtile.sampleType = sampleType;

  if (loglevel >= 5)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  END" << endl;
  return rawtile;
}


// Single display settings wrapper around general display settings.
class SingleDisplaySettingsWrapper : public libCZI::IDisplaySettings
{
  // Wrapper around general channel display setting, with GetIsEnabled() forced to be true.
  class EnabledChannelDisplaySetting : public libCZI::IChannelDisplaySetting
  {
  private:
    const std::shared_ptr<libCZI::IChannelDisplaySetting> enabled_channel_display_setting;
  public:
    explicit EnabledChannelDisplaySetting(const std::shared_ptr<libCZI::IChannelDisplaySetting>& channel_display_setting)
      : enabled_channel_display_setting(channel_display_setting)
    {}

    std::shared_ptr<libCZI::IChannelDisplaySetting> GetChannelDisplaySetting() const {
      return enabled_channel_display_setting;
    }

    // Override GetIsEnabled() to always return true.
    bool GetIsEnabled() const override {
      return true;
    }

    // Pass through remaining libCZI::IChannelDisplaySetting methods.
    float    GetWeight() const override {
      return enabled_channel_display_setting->GetWeight();
    }
    bool    TryGetTintingColorRgb8(libCZI::Rgb8Color* pColor) const override {
      return enabled_channel_display_setting->TryGetTintingColorRgb8(pColor);
    }
    void    GetBlackWhitePoint(float* pBlack, float* pWhite) const override {
      enabled_channel_display_setting->GetBlackWhitePoint(pBlack, pWhite);
    }
    libCZI::IDisplaySettings::GradationCurveMode GetGradationCurveMode() const override {
      return enabled_channel_display_setting->GetGradationCurveMode();
    }
    bool    TryGetGamma(float* gamma) const override {
      return enabled_channel_display_setting->TryGetGamma(gamma);
    }
    bool    TryGetSplineControlPoints(std::vector<libCZI::IDisplaySettings::SplineControlPoint>* ctrlPts) const override {
      return enabled_channel_display_setting->TryGetSplineControlPoints(ctrlPts);
    }
    bool    TryGetSplineData(std::vector<libCZI::IDisplaySettings::SplineData>* data) const override {
      return enabled_channel_display_setting->TryGetSplineData(data);
    }
  };

private:
  const EnabledChannelDisplaySetting wrapped_display_setting;
public:
  explicit SingleDisplaySettingsWrapper(const std::shared_ptr<libCZI::IDisplaySettings>& display_settings,
                                        int channel_index)
    : wrapped_display_setting(display_settings->GetChannelDisplaySettings(channel_index))
  {}

  void EnumChannels(std::function<bool(int chIndex)> func) const override {
    func(0);
  }

  std::shared_ptr<libCZI::IChannelDisplaySetting> GetChannelDisplaySettings(int chIndex) const override {
    if (chIndex == 0)
      return std::make_shared<EnabledChannelDisplaySetting>(wrapped_display_setting);

    return std::shared_ptr<libCZI::IChannelDisplaySetting>();
  }
};


// See:
//    MultiChannelCompositor.cpp -- class CMultiChannelCompositor2
//        -- static void ComposeMultiChannel_Bgr24()
//        -- static void ComposeMultiChannel()
//        -- struct FunctionsBgr24 {}
//        -- static void DoTintingBlackWhitePtCopy()
//        -- static void DoTintingBlackWhitePt()
//        -- static void CopyTinting()
//        -- struct CGetBlackWhitePtTintingGray16 {}
//        -- struct CGetBlackWhitePtGray16 {}

struct PixelGraderBase {
  const libCZI::Compositors::ChannelInfo channelInfo;

  PixelGraderBase(const libCZI::Compositors::ChannelInfo& channel_info)
    : channelInfo(channel_info) {}

  virtual void CopyPt(uint8_t* dst, const uint8_t* src) const = 0;
};

template <typename tValue, int maxValue>
struct PixelGrader : PixelGraderBase {
  PixelGrader(const libCZI::Compositors::ChannelInfo& channel_info)
    : PixelGraderBase(channel_info) {
    if (loglevel >= 9)
      logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
              << ", maxValue " << maxValue
              << ", blackPoint " << channelInfo.blackPoint << " > 0.0 ? " << (channelInfo.blackPoint > 0.0)
              << ", whitePoint " << channelInfo.whitePoint << " < 1.0 ? " << (channelInfo.whitePoint < 1.0)
              << endl;
  }

  void CopyPt(uint8_t* dst, const uint8_t* src) const {
    float src_float = (float) *((tValue *) src) / (float) maxValue;
    float dst_float;
    if (src_float <= channelInfo.blackPoint) {
      dst_float = 0.0;
    }
    else if (src_float >= channelInfo.whitePoint) {
      dst_float = 1.0;
    }
    else {
      dst_float = (src_float - channelInfo.blackPoint) /
        (channelInfo.whitePoint - channelInfo.blackPoint);
    }

    *((tValue *) dst) = (tValue)(dst_float * maxValue + .5);
  }
};

shared_ptr<libCZI::IBitmapData> GradeSingleChannel_Gray(const shared_ptr<libCZI::IBitmapData>& src_bitmap,
                                                        const libCZI::Compositors::ChannelInfo& channelInfo,
                                                        const libCZI::RgbFloatColor& background) {
  shared_ptr<libCZI::IBitmapData> dst_bitmap =
    GetSite()->CreateBitmap(src_bitmap->GetPixelType(), src_bitmap->GetWidth(), src_bitmap->GetHeight());

  libCZI::ScopedBitmapLockerP locked_src{src_bitmap.get()};
  libCZI::ScopedBitmapLockerP locked_dst{dst_bitmap.get()};

  shared_ptr<PixelGraderBase> pixel_grader;
  if (src_bitmap->GetPixelType() == libCZI::PixelType::Gray16
      && (0.0 < channelInfo.blackPoint || channelInfo.whitePoint < 1.0)) {
    pixel_grader = shared_ptr<PixelGraderBase>( new PixelGrader<std::uint16_t, 256 * 256 - 1>( channelInfo ) );
  }
  else if (src_bitmap->GetPixelType() == libCZI::PixelType::Gray8
           && (0.0 < channelInfo.blackPoint || channelInfo.whitePoint < 1.0)) {
    pixel_grader = shared_ptr<PixelGraderBase>( new PixelGrader<std::uint8_t, 256 - 1>( channelInfo ) );
  }
  else {
    // Just copy without any gradation.
    if (loglevel >= 9)
      logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  Copy w/o gradation" << endl;

    CBitmapOperations::Fill(dst_bitmap.get(), background);  // In case of non-copied default background.

    CBitmapOperations::Copy(src_bitmap->GetPixelType(), locked_src.ptrDataRoi, locked_src.stride,
                            dst_bitmap->GetPixelType(), locked_dst.ptrDataRoi, locked_dst.stride,
                            src_bitmap->GetWidth(), src_bitmap->GetHeight(), false);

    return dst_bitmap;
  }


  uint32_t w = src_bitmap->GetWidth();
  uint32_t h = src_bitmap->GetHeight();
  uint8_t bytes_per_pel = CziUtils::GetBytesPerPel(src_bitmap->GetPixelType());

  const uint8_t* ptrSrc = (uint8_t*) locked_src.ptrDataRoi;
  uint32_t strideSrc= locked_src.stride;

  uint8_t* ptrDst = (uint8_t*) locked_dst.ptrDataRoi;
  uint32_t strideDst= locked_dst.stride;

  for (uint32_t y = 0; y < h; ++y) {
    const uint8_t* pSrc = ptrSrc + y * strideSrc;
    uint8_t* pDst = ptrDst + y * strideDst;

    for (uint32_t x = 0; x < w; ++x) {
      pixel_grader->CopyPt(pDst, pSrc);

      pSrc += bytes_per_pel;
      pDst += bytes_per_pel;
    }
  }

  return dst_bitmap;
}


// Multiple CZI channels (usu. fluorescence, Gray16).
RawTile CZIImage::getAllChannelsPyramidLayerTile(
    const int seq, const int ang, const unsigned int res, const unsigned int tile,
    const int z_layer, const int czi_pyr_layer,
    const unsigned int tile_w, const unsigned int tile_h, const libCZI::IntRect roi)
{
  if (loglevel >= 2)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__
            << "( seq " << seq << ", ang " << ang << ", res " << res << ", tile " << tile
            << ", z_layer " << z_layer << ", czi_pyr_layer " << czi_pyr_layer
            << ", tile_w " << tile_w << ", tile_h " << tile_h << ", roi " << roi
            << " )  BEGIN" << endl;

  /// For fluorescence images, using Command::ScalingChannelComposite works well.
  // See:  libCZI/Src/CZICmd/execute.cpp -- class CExecuteScalingChannelComposite -- execute().

  libCZI::CDimCoordinate coordinate;
  if (z_layers_size > 0) {
    coordinate.Set(libCZI::DimensionIndex::Z, z_layers_start + z_layer);
  }

  libCZI::ISingleChannelPyramidLayerTileAccessor::PyramidLayerInfo pyr_layer_info;
  pyr_layer_info.minificationFactor = image_minification;
  pyr_layer_info.pyramidLayerNo = czi_pyr_layer;

  libCZI::ISingleChannelPyramidLayerTileAccessor::Options scpta_options;
  scpta_options.Clear();
  //  scpta_options.sceneFilter = options.GetSceneIndexSet(); // Unused, leave as default from Clear().
  libCZI::RgbFloatColor bright_bkgd{ 0.9, 0.9, 0.9 };  // Light for brightfield images.
  libCZI::RgbFloatColor fluor_bkgd{ 0.0, 0.0, 0.0 };  // Black for fluorescence channels.
  scpta_options.backGroundColor = (channels_size == 1 ? bright_bkgd : fluor_bkgd);
  /*TEMP(for subblock debugging)  scpta_options.drawTileBorder = true;*/


  std::shared_ptr<libCZI::IDisplaySettings> full_display_settings = (czi_reader
                                                                     ->ReadMetadataSegment()
                                                                     ->CreateMetaFromMetadataSegment()
                                                                     ->GetDocumentInfo()
                                                                     ->GetDisplaySettings());
  // Note: IIPImage::channels, the base class number of channels,
  // must be set in CZIImage::loadImageInfo(),
  // and in either #ifdef case below,
  // we must end up with:  channels == active_channels.size()
#ifdef USE_ONLY_ACTIVE_ENABLED_CHANNELS
  std::vector<int> active_channels = libCZI::CDisplaySettingsHelper::GetActiveChannels(full_display_settings.get());
#else  // Use all channels
  std::vector<int> active_channels;
  for (int idx = 0; idx < channels; ++idx) { active_channels.push_back(idx); }
#endif

  auto accessor = czi_reader->CreateSingleChannelPyramidLayerTileAccessor();
  std::vector<shared_ptr<libCZI::IBitmapData>> channel_bitmaps;

  for (int idx = 0; idx < (int) channels; ++idx) {
    int channel_num = active_channels.at(idx);

    if (czi_reader->GetStatistics().dimBounds.IsValid(libCZI::DimensionIndex::C)) {
      // That's a cornerstone case - or a loophole in the specification: if the document
      // does not contain C-dimension (=none of the sub-blocks has a valid C-dimension),
      // then we must not set the C-dimension here. I suppose we should define that a
      // valid C-dimension is mandatory...
      coordinate.Set(libCZI::DimensionIndex::C, channel_num);
    }

    channel_bitmaps.emplace_back(accessor->Get(roi, &coordinate, pyr_layer_info, &scpta_options));
  }


  std::uint32_t tile_buf_width = channel_bitmaps[0]->GetWidth();
  std::uint32_t tile_buf_height = channel_bitmaps[0]->GetHeight();
  tile_buf_malloc(tile_buf_width * tile_buf_height * channels * (bpc/8));


  if (loglevel >= 2)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << ", tile_buf_size " << tile_buf_size
            << ", tile_buf_width " << tile_buf_width << ", tile_buf_height " << tile_buf_height
            << ", channels_size " << channels_size
            << ", active_channels.size() " << active_channels.size()
            << ", bpc " << bpc << ", bpc/8 " << bpc/8
            << ", tile_w " << tile_w << ", tile_h " << tile_h
            << endl;


  // Copy bitmaps channel-by-channel, collating them into a single channels-wide bitmap.
  for (int idx = 0; idx < (int) channels; ++idx) {
    int channel_num = active_channels.at(idx);

    std::shared_ptr<libCZI::IDisplaySettings> single_display_settings(
        new SingleDisplaySettingsWrapper(full_display_settings, channel_num));

    libCZI::CDisplaySettingsHelper single_settings_helper;
    single_settings_helper.Initialize(single_display_settings.get(),
                                      [&](int chIndx)->libCZI::PixelType {
                                        // Note: idx is the only channel index used.
                                        return channel_bitmaps[idx]->GetPixelType();
                                      });

    libCZI::Compositors::ChannelInfo single_channel_info =
      single_settings_helper.GetChannelInfosArray()[0];

    if (loglevel >= 2)
      logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
              << ", idx " << idx << ", channel " << channel_num
              << ": weight " << single_channel_info.weight
              << ", original enableTinting " << single_channel_info.enableTinting
              << ", tinting bgr " << (int)single_channel_info.tinting.color.b
              << "." << (int)single_channel_info.tinting.color.g
              << "." << (int)single_channel_info.tinting.color.r
              << ", blackPoint " << single_channel_info.blackPoint
              << ", whitePoint " << single_channel_info.whitePoint
              << endl;

    // Manually alter channel info to full gray-scale (#FFFFFF) tinting
    // by disabling tinting, rather than tinting with { 255, 255, 255 }.
    single_channel_info.enableTinting = false;

    shared_ptr<libCZI::IBitmapData> graded_bitmap = GradeSingleChannel_Gray(channel_bitmaps[idx],
                                                                            single_channel_info,
                                                                            scpta_options.backGroundColor);
    libCZI::ScopedBitmapLockerP locked_graded{graded_bitmap.get()};

    if (loglevel >= 2)
      logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
              << ", locked_graded.size " << locked_graded.size
              << ", locked_graded.stride " << locked_graded.stride
              << ", channel_bitmaps[" << idx << "]->GetWidth() " << channel_bitmaps[idx]->GetWidth()
              << " <=> graded_bitmap->GetWidth() " << graded_bitmap->GetWidth()
              << ", channel_bitmaps[" << idx << "]->GetHeight() " << channel_bitmaps[idx]->GetHeight()
              << " <=> graded_bitmap->GetHeight() " << graded_bitmap->GetHeight()
              << endl;


    // Collate graded channel into channels == active_channels.size() deep tile buffer.
    std::uint32_t graded_px_width = locked_graded.stride / (bpc/8);  // Stride is in bytes.
    for (std::uint32_t h = 0; h < tile_buf_height; ++h) {
      for (std::uint32_t w = 0; w < tile_buf_width; ++w) {
        std::uint32_t tile_px = h * tile_buf_width + w;
        std::uint32_t graded_px = h * graded_px_width + w;

        if (bpc == 32) {
          if (sampleType == FLOATINGPOINT) {
            ((float *) tile_buf)[tile_px * channels + idx] =
              ((float *) locked_graded.ptrDataRoi)[graded_px];
          }
          else {
            ((std::uint32_t *) tile_buf)[tile_px * channels + idx] =
              ((std::uint32_t *) locked_graded.ptrDataRoi)[graded_px];
          }
        }
        else if (bpc == 16) {
          ((std::uint16_t *) tile_buf)[tile_px * channels + idx] =
            ((std::uint16_t *) locked_graded.ptrDataRoi)[graded_px];
        }
        else {
          ((std::uint8_t *) tile_buf)[tile_px * channels + idx] =
            ((std::uint8_t *) locked_graded.ptrDataRoi)[graded_px];
        }
      }
    }
  }


  if (loglevel >= 2)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()"
            << "  RawTile(): tile " << tile << ", res " << res
            << ", seq " << seq << ", ang " << ang
            << ", tile_w " << tile_w << ", tile_h " << tile_h
            << ", channels " << channels << ", bpc " << bpc
            << endl;

  RawTile rawtile( tile, res, seq, ang, tile_w, tile_h, channels, bpc );
  rawtile.data = tile_buf;
  rawtile.dataLength = tile_buf_size;
  rawtile.filename = getImagePath();
  rawtile.timestamp = timestamp;
  rawtile.memoryManaged = 0;
  //rawtile.padded = true;  // TODO(Leo) Huh? Why padded for QPTIFF, not for OpenSlide?
  rawtile.sampleType = sampleType;

  if (loglevel >= 5)
    logfile << __FILE__ << ": " << __LINE__ << "  " << __FUNCTION__ << "()  END" << endl;
  return rawtile;
}
