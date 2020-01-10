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

  // Try to open and allocate a CZI-reader.
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



void CZIImage::loadImageInfo( int x, int y )
{
  currentX = x;
  currentY = y;

  auto mds = czi_reader->ReadMetadataSegment();
  auto md = mds->CreateMetaFromMetadataSegment();

  auto sbStatistics = czi_reader->GetStatistics();

  // Store the list of image dimensions available
  // See:  CZICmd/execute.cpp -- PrintStatistics().
  ////unsigned int full_width = sbStatistics.boundingBox.w;
  ////unsigned int full_height = sbStatistics.boundingBox.h;
  unsigned int full_width = sbStatistics.boundingBoxLayer0Only.w;
  unsigned int full_height = sbStatistics.boundingBoxLayer0Only.h;
  image_widths.push_back( full_width );
  image_heights.push_back( full_height );
  numResolutions = 1;

  // CZI pyramid layers as virtual width x height images.
  // See:  CZICmd/execute.cpp -- PrintPyramidStatistics().
  auto pyrStatistics = czi_reader->GetPyramidStatistics();
  auto scene0PyrStatistics = pyrStatistics.scenePyramidStatistics[0];
  for (const auto& j : scene0PyrStatistics) {
    if (!j.layerInfo.IsNotIdentifiedAsPyramidLayer()) {
      if (j.layerInfo.IsLayer0() != true) {
        int scale = j.layerInfo.minificationFactor;
        for (int n = 1; n < j.layerInfo.pyramidLayerNo; ++n) {
          scale *= j.layerInfo.minificationFactor;
        }

		// Scaled down dimensions [rounded up with (.. -1)/scale +1]
		unsigned int w = (full_width -1)/scale +1;
		unsigned int h = (full_height -1)/scale +1;
		// Ignore downsamples smaller than 2K x 2K.  TODO(Leo) Why?  Ask Coleman.
		if (w < 2000 && h < 2000) {
		  break;
		}

		image_widths.push_back( w );
		image_heights.push_back( h );
		++numResolutions;
      }
    }
  }
#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// The image pixel dimensions
  std::vector <unsigned int> image_widths, image_heights;
  /*
	static void PrintStatistics(const CCmdLineOptions& options, ICZIReader* reader)
		WriteIntRect(ss, sbStatistics.boundingBox);


18_WT1_NKPP002_01_R_LS.czi.PrintInformation.txt
Bounding-Box:
 All:    X=0 Y=0 W=98172 H=66096
 Layer0: X=0 Y=0 W=98097 H=65963

scene#0:
 number of subblocks with scale 1/1: 2550
 number of subblocks with scale 1/3: 919
 number of subblocks with scale 1/9: 116
 number of subblocks with scale 1/27: 17
 number of subblocks with scale 1/81: 3
 number of subblocks with scale 1/243: 1

>>> 98172/243.
404.0
>>> 66096/243.
272.0


GR57-13 2015_10_12__0078.czi.PrintInformation.txt
Bounding-Box:
 All:    X=-134078 Y=9180 W=119630 H=48000
 Layer0: X=-134078 Y=9180 W=119625 H=47998

scene#0:
 number of subblocks with scale 1/1: 1311
 number of subblocks with scale 1/2: 1221
 number of subblocks with scale 1/4: 333
 number of subblocks with scale 1/8: 99
 number of subblocks with scale 1/16: 30
 number of subblocks with scale 1/32: 12
 number of subblocks with scale 1/64: 3

scene#1:
 number of subblocks with scale 1/1: 1182
 number of subblocks with spcale 1/2: 1110
 number of subblocks with scale 1/4: 303
 number of subblocks with scale 1/8: 90
 number of subblocks with scale 1/16: 33
 number of subblocks with scale 1/32: 12
 number of subblocks with scale 1/64: 3

>>> 119630/64
1869.21875
>>> 48000/64
750.0
>>> 119630/4 
29907.5


1_CD3645_Ins594_Glu488_DAPI_IS4612.czi.PrintInformation.txt
Bounding-Box:
 All:    X=-111996 Y=34869 W=5731 H=7580
 Layer0: X=-111996 Y=34869 W=5731 H=7580

scene#0:
 number of subblocks with scale 1/1: 480
 number of subblocks with scale 1/2: 480
 number of subblocks with scale 1/4: 160
 number of subblocks with scale 1/8: 40

>>> 7580/8
947.5
>>> 5731/8
716.375
>>> 7580/4
1895.0
>>> 5731/4
1432.75
  */
#endif // TODO(Leo)  IIPImage metadata to be initialized.

  // CZI images are not stored as grids of tiles;
  // virtual tiles will be created on the fly.
  tile_width = 1024;
  tile_height = 1024;
#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// The base tile pixel dimensions
  unsigned int tile_width, tile_height;
  /*
CZI aren't arranged as tiles, so output a default virtual tile dimensions, e.g.: 728x728 or 1024x1024
  */

#endif // TODO(Leo)  IIPImage metadata to be initialized.

#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// The colour space of the image
  ColourSpaces colourspace;
  /*
execute.cpp:: static void PrintDisplaySettingsForChannel(int ch, IChannelDisplaySetting* dsplChSettings, const CCmdLineOptions& options)

Complete list of sub-blocks
---------------------------

#0: C0S0B0 M=0 logical=(0,0,98172,66096) phys.=(404,272) pixeltype=bgr24 comp.mode=jpgxr


$ grep -Po 'pixeltype=\S+|^ C -> Start=\d+ Size=\d+' -r | sort | uniq -c
      1 18_WT1_NKPP002_01_R_LS.czi.PrintInformation.txt: C -> Start=0 Size=1
   3606 18_WT1_NKPP002_01_R_LS.czi.PrintInformation.txt:pixeltype=bgr24
      1 1_CD3645_Ins594_Glu488_DAPI_IS4612.czi.PrintInformation.txt: C -> Start=0 Size=4
   1160 1_CD3645_Ins594_Glu488_DAPI_IS4612.czi.PrintInformation.txt:pixeltype=gray16
  */


#endif // TODO(Leo)  IIPImage metadata to be initialized.

#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// The number of available resolutions in this image
  // CZIcmd --source ../AC5537\ \ insulin\ CD4.czi --command PrintInformation --info-level 'PyramidStatistics'
  unsigned int numResolutions;  // --command PrintInformation --info-level 'PyramidStatistics'
  /*
	static void PrintPyramidStatistics(ICZIReader* reader, const CCmdLineOptions& options)

      1 18_WT1_NKPP002_01_R_LS.czi.PrintInformation.txt: C -> Start=0 Size=1
  */

#endif // TODO(Leo)  IIPImage metadata to be initialized.

#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// The bits per channel for this image
  unsigned int bpc;

#endif // TODO(Leo)  IIPImage metadata to be initialized.

#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// The number of channels for this image
  // CZIcmd --source ../1_CD3645_Ins594_Glu488_DAPI_IS4612.czi --command PrintInformation --info-level 'Statistics'
  //     C -> Start=0 Size=4
  unsigned int channels;

#endif // TODO(Leo)  IIPImage metadata to be initialized.

#if 0  // TODO(Leo)  IIPImage metadata to be initialized.
  /// The sample format type (fixed or floating point)
  SampleType sampleType;

#endif // TODO(Leo)  IIPImage metadata to be initialized.

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

  /// Indicate whether we have opened and initialised some parameters for this image
  //  bool isSet;  // Done - CZIImage::openImage()

  /// If we have an image sequence, the current X and Y position
  //  int currentX, currentY;  // Below

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

  /// Image modification timestamp
  //  time_t timestamp;  // Done - CZIImage::openImage() .. updateTimestamp()
#endif // TODO(Leo)

#if 0  // TODO(Leo)
  tdir_t current_dir;
  int count;
  uint16 colour, samplesperpixel, bitspersample, sampleformat;
  double sminvaluearr[4] = {0.0}, smaxvaluearr[4] = {0.0};
  double *sminvalue = NULL, *smaxvalue = NULL;
  unsigned int w, h;
  string filename;
  char *tmp = NULL;

  currentX = x;
  currentY = y;

  // Get the tile and image sizes
#if 0  // TODO(Leo)
  // TODO(Leo)
  /*
	CZICmd for appropriate metadata???
	    Extract width and height from CZIcmd Statistics!!!
CZIcmd --source "${DOWNLOAD}" \
		   --command PrintInformation --info-level 'Statistics'
  */
#endif // TODO(Leo)
  CZIGetField( czi, CZITAG_TILEWIDTH, &tile_width );
  CZIGetField( czi, CZITAG_TILELENGTH, &tile_height );
  CZIGetField( czi, CZITAG_IMAGEWIDTH, &w );
  CZIGetField( czi, CZITAG_IMAGELENGTH, &h );
  CZIGetField( czi, CZITAG_SAMPLESPERPIXEL, &samplesperpixel );
  CZIGetField( czi, CZITAG_BITSPERSAMPLE, &bitspersample );
  CZIGetField( czi, CZITAG_PHOTOMETRIC, &colour );
  CZIGetField( czi, CZITAG_SAMPLEFORMAT, &sampleformat );

  if ( samplesperpixel == 1 ){
    channels = 0;

    unsigned int w1 = w, h1 = h;
    while ( w1 == w && h1 == h ){
      CZIReadDirectory( czi );
      CZIGetField( czi, CZITAG_IMAGEWIDTH, &w1 );
      CZIGetField( czi, CZITAG_IMAGELENGTH, &h1 );
      channels++;
      logfile << "CZIImage:: " << "Dir w " << w1 << " Chan " << channels << endl;
    }
  }
  else {
    channels = 3;
  }

  logfile << "CZIImage:: " << "Channels: " << channels << endl;

  bpc = (unsigned int) bitspersample;
  sampleType = (sampleformat==3) ? FLOATINGPOINT : FIXEDPOINT;

  // // Check for the no. of resolutions in the pyramidal image
  // current_dir = CZICurrentDirectory( czi );
  // CZISetDirectory( czi, 0 );

  // Store the list of image dimensions available
  image_widths.push_back( w );
  image_heights.push_back( h );

  for( count = 0; CZIReadDirectory( czi ); count++ ){
    CZIGetField( czi, CZITAG_IMAGEWIDTH, &w );
    CZIGetField( czi, CZITAG_IMAGELENGTH, &h );
    // ignore downsamples smaller than 2K x 2K
    if (w < 2000 && h < 2000) {
      break;
    }
    image_widths.push_back( w );
    image_heights.push_back( h );
    current_dir = CZICurrentDirectory( czi );
    CZISetDirectory( czi, current_dir + channels );
    logfile << "CZIImage:: " << "Channels: " << count << " w " << w << endl;
  }
  // Reset the CZI directory
  CZISetDirectory( czi, current_dir );

  numResolutions = count+1;

  // Handle various colour spaces
  if( colour == PHOTOMETRIC_CIELAB ) colourspace = CIELAB;
  else if( colour == PHOTOMETRIC_MINISBLACK ){
    colourspace = (bpc==1)? BINARY : GREYSCALE;
  }
  else if( colour == PHOTOMETRIC_PALETTE ){
    // Watch out for colourmapped images. These are stored as 1 sample per pixel,
    // but are decoded to 3 channels by libczi, so declare them as sRGB
    colourspace = sRGB;
    channels = 3;
  }
  else if( colour == PHOTOMETRIC_YCBCR ){
    // JPEG encoded tiles can be subsampled YCbCr encoded. Ask to decode these to RGB
    CZISetField( czi, CZITAG_JPEGCOLORMODE, JPEGCOLORMODE_RGB );
    colourspace = sRGB;
  }
  else colourspace = sRGB;

  // Get the max and min values for our data type - required for floats
  // This are usually single values per image, but can also be per channel
  // in libczi > 4.0.2 via http://www.asmail.be/msg0055458208.html

#ifdef CZITAG_PERSAMPLE
  if( channels > 1 ){
    CZISetField(czi, CZITAG_PERSAMPLE, PERSAMPLE_MULTI);
    CZIGetFieldDefaulted( czi, CZITAG_SMINSAMPLEVALUE, &sminvalue );
    CZIGetFieldDefaulted( czi, CZITAG_SMAXSAMPLEVALUE, &smaxvalue );
    CZISetField(czi, CZITAG_PERSAMPLE, PERSAMPLE_MERGED);
    if (!sminvalue) sminvalue = sminvaluearr;
    if (!smaxvalue) smaxvalue = smaxvaluearr;
  }
  else{
#endif
    sminvalue = sminvaluearr;
    smaxvalue = smaxvaluearr;
    CZIGetFieldDefaulted( czi, CZITAG_SMINSAMPLEVALUE, sminvalue );
    CZIGetFieldDefaulted( czi, CZITAG_SMAXSAMPLEVALUE, smaxvalue );
#ifdef CZITAG_PERSAMPLE
  }
#endif

  // Clear our arrays
  min.clear();
  max.clear();

  for( unsigned int i=0; i<channels; i++ ){
    if( (!sminvalue) == smaxvalue[i] ){
      // Set default values if values not included in header
      if( bpc <= 8 ) smaxvalue[i] = 255.0;
      else if( bpc == 12 ) smaxvalue[i] = 4095.0;
      else if( bpc == 16 ) smaxvalue[i] = 65535.0;
      else if( bpc == 32 && sampleType == FIXEDPOINT ) smaxvalue[i] = 4294967295.0;
    }
    min.push_back( (float)sminvalue[i] );
    max.push_back( (float)smaxvalue[i] );
  }

  // Also get some basic metadata
  if( CZIGetField( czi, CZITAG_ARTIST, &tmp ) ) metadata["author"] = tmp;
  if( CZIGetField( czi, CZITAG_COPYRIGHT, &tmp ) ) metadata["copyright"] = tmp;
  if( CZIGetField( czi, CZITAG_DATETIME, &tmp ) ) metadata["create-dtm"] = tmp;
  if( CZIGetField( czi, CZITAG_IMAGEDESCRIPTION, &tmp ) ) metadata["subject"] = tmp;
  if( CZIGetField( czi, CZITAG_SOFTWARE, &tmp ) ) metadata["app-name"] = tmp;
  if( CZIGetField( czi, CZITAG_XMLPACKET, &count, &tmp ) ) metadata["xmp"] = string(tmp,count);
  if( CZIGetField( czi, CZITAG_ICCPROFILE, &count, &tmp ) ) metadata["icc"] = string(tmp,count);

#endif // TODO(Leo)
}

RawTile CZIImage::getTile( int x, int y, unsigned int res, int layers, unsigned int tile )
{
  return RawTile();

#if 0  // TODO(Leo)

  uint32 im_width, im_height, tw, th, ntlx, ntly;
  uint32 rem_x, rem_y;
  uint16 colour;
  string filename;


  logfile << "CZIIMAGE:: getTile begin" << endl;

  // Check the resolution exists
  if( res > numResolutions ){
    ostringstream error;
    error << "CZIImage :: Asked for non-existent resolution: " << res;
    throw file_error( error.str() );
  }


  // If we are currently working on a different sequence number, then
  //  close and reload the image.
  if( (currentX != x) || (currentY != y) ){
    closeImage();
  }


  // Open the CZI if it's not already open
  if( !czi ){
    filename = getFileName( x, y );
    if( ( czi = CZIOpen( filename.c_str(), "rm" ) ) == NULL ){
      throw file_error( "czi open failed for:" + filename );
    }
  }


  // Reload our image information in case the tile size etc is different
  if( (currentX != x) || (currentY != y) ){
    loadImageInfo( x, y );
  }


  // The first resolution is the highest, so we need to invert
  //  the resolution - can avoid this if we store our images with
  //  the smallest image first.
  int vipsres = ( numResolutions - 1 ) - res;

  // because we have N directories per resolution and because the N+1th directory is the thumbnail
  int czi_dir = vipsres * channels;
  if ( czi_dir > 0 ) {
    czi_dir += 1;
  }

  logfile << "CZIIMAGE:: czi_dir = " << czi_dir << " bpc = " << bpc << endl;

  // Change to the right directory for the resolution
  if( !CZISetDirectory( czi, czi_dir ) ) {
    throw file_error( "CZISetDirectory failed" );
  }


  // Check that a valid tile number was given
  if( tile >= CZINumberOfTiles( czi ) ) {
    ostringstream tile_no;
    tile_no << "Asked for non-existent tile: " << tile;
    throw file_error( tile_no.str() );
  }


  // Get the size of this tile, the current image,
  //  the number of samples and the colourspace.
  // CZITAG_TILEWIDTH give us the values for the resolution,
  //  not for the tile itself
  CZIGetField( czi, CZITAG_TILEWIDTH, &tw );
  CZIGetField( czi, CZITAG_TILELENGTH, &th );
  CZIGetField( czi, CZITAG_IMAGEWIDTH, &im_width );
  CZIGetField( czi, CZITAG_IMAGELENGTH, &im_height );
  CZIGetField( czi, CZITAG_PHOTOMETRIC, &colour );
//   CZIGetField( czi, CZITAG_SAMPLESPERPIXEL, &channels );
//   CZIGetField( czi, CZITAG_BITSPERSAMPLE, &bpc );


  // Make sure this resolution is tiled
  if( (tw == 0) || (th == 0) ){
    throw file_error( "Requested resolution is not tiled" );
  }


  // Total number of bytes in tile
  unsigned int np = tw * th;


  // Get the width and height for last row and column tiles
  rem_x = im_width % tw;
  rem_y = im_height % th;


  // Calculate the number of tiles in each direction
  ntlx = (im_width / tw) + (rem_x == 0 ? 0 : 1);
  ntly = (im_height / th) + (rem_y == 0 ? 0 : 1);


  // Alter the tile size if it's in the last column
  if( ( tile % ntlx == ntlx - 1 ) && ( rem_x != 0 ) ) {
    tw = rem_x;
  }


  // Alter the tile size if it's in the bottom row
  if( ( tile / ntlx == ntly - 1 ) && rem_y != 0 ) {
    th = rem_y;
  }


  // Handle various colour spaces
  if( colour == PHOTOMETRIC_CIELAB ) colourspace = CIELAB;
  else if( colour == PHOTOMETRIC_MINISBLACK ){
    colourspace = (bpc==1)? BINARY : GREYSCALE;
  }
  else if( colour == PHOTOMETRIC_PALETTE ){
    // Watch out for colourmapped images. There are stored as 1 sample per pixel,
    // but are decoded to 3 channels by libczi, so declare them as sRGB
    colourspace = GREYSCALE;
    channels = 1;
  }
  else if( colour == PHOTOMETRIC_YCBCR ){
    // JPEG encoded tiles can be subsampled YCbCr encoded. Ask to decode these to RGB
    CZISetField( czi, CZITAG_JPEGCOLORMODE, JPEGCOLORMODE_RGB );
    colourspace = sRGB;
  }
  else colourspace = sRGB;

  logfile << "CZIImage:: CZITileSize = " << CZITileSize(czi) << endl;
  tsize_t tile_size = CZITileSize(czi);

  tdata_t *channels_tile_bufs = new tdata_t[channels];
  for ( int i = 0; i < channels; i++ ) {
    if( ( channels_tile_bufs[i] = _CZImalloc ( tile_size ) ) == NULL ) {
      throw file_error( "czi malloc tile failed" );
    }
    int length = CZIReadEncodedTile( czi, (ttile_t) tile,
            channels_tile_bufs[i], (tsize_t) - 1 );
    if( length == -1 ) {
      throw file_error( "CZIReadEncodedTile failed for " + getFileName( x, y ) );
    }
    CZIReadDirectory( czi );
  }   

  // Allocate memory for our tile.
  if( !tile_buf ){
    logfile << "CZIIMAGE:: tile_buf size = " << tile_size * channels << endl;
    if( ( tile_buf = _CZImalloc( tile_size * channels ) ) == NULL ){
      throw file_error( "czi malloc tile failed" );
    }
  }

  logfile << "CZIIMAGE:: tiles read into channels_tile_bufs" << endl;

  // Decode and read the tile
  // int length = CZIReadEncodedTile( czi, (ttile_t) tile,
  //   tile_buf, (tsize_t) - 1 );
  // if( length == -1 ) {
  //   throw file_error( "CZIReadEncodedTile failed for " + getFileName( x, y ) );
  // }

  for ( int i = 0; i < np; i++ ){
    for ( int c = 0; c < channels; c++ ){
      if (bpc == 32) {
        if (sampleType == FLOATINGPOINT) {
          ((float *)tile_buf)[i * channels + c] = ((float *)channels_tile_bufs[c])[i];
        }
        else {
          ((unsigned int *)tile_buf)[i * channels + c] = ((unsigned int *)channels_tile_bufs[c])[i];
        }
      }
      else if (bpc == 16) {
        ((unsigned short *)tile_buf)[i * channels + c] = ((unsigned short *)channels_tile_bufs[c])[i];
      }
      else {
        ((unsigned char *)tile_buf)[i * channels + c] = ((unsigned char *)channels_tile_bufs[c])[i];
      }
    }
  }

  for ( int i = 0; i < channels; i++ ){
    _CZIfree(channels_tile_bufs[i]);
  }
  delete channels_tile_bufs;

  logfile << " CZIImage:: bpc = " << bpc << endl;

  RawTile rawtile( tile, res, x, y, tw, th, channels, bpc );
  rawtile.data = tile_buf;
  rawtile.dataLength = tile_size * channels;
  rawtile.filename = getImagePath();
  rawtile.timestamp = timestamp;
  rawtile.memoryManaged = 0;
  rawtile.padded = true;
  rawtile.sampleType = sampleType;


  // Pad 1 bit 1 channel bilevel images to 8 bits for output
  if( bpc==1 && channels==1 ){

    // Pixel index
    unsigned int n = 0;

    // Calculate number of bytes used - round integer up efficiently
    unsigned int nbytes = (np + 7) / 8;
    unsigned char *buffer = new unsigned char[np];

    // Take into account photometric interpretation:
    //   0: white is zero, 1: black is zero
    unsigned char min = (unsigned char) 0;
    unsigned char max = (unsigned char) 255;
    if( colour == 0 ){
      min = (unsigned char) 255; max = (unsigned char) 0;
    }

    // Unpack each raw byte into 8 8-bit pixels
    for( unsigned int i=0; i<nbytes; i++ ){
      unsigned char t = ((unsigned char*)tile_buf)[i];
      // Count backwards as CZI is usually MSB2LSB
      for( int k=7; k>=0; k-- ){
        // Set values depending on whether bit is set
        buffer[n++] = (t & (1 << k)) ? max : min;
      }
    }

    rawtile.dataLength = n;
    rawtile.data = buffer;
    rawtile.bpc = 8;
    rawtile.memoryManaged = 1;
  }


  return( rawtile );

#endif // TODO(Leo)

}
