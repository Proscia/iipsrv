// Member functions for QPTIFFImage.h

/*  IIP Server: Tiled Pyramidal TIFF handler

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


#include "QPTIFFImage.h"
#include <sstream>
#include <iostream>
#include <fstream>

extern int loglevel;
extern std::ofstream logfile;

using namespace std;


void QPTIFFImage::openImage()
{

  // Insist that the tiff and tile_buf be NULL
  if( tiff || tile_buf ){
    throw file_error( "TPT::openImage: tiff or tile_buf is not NULL" );
  }

  string filename = getFileName( currentX, currentY );

  // Update our timestamp
  updateTimestamp( filename );

  // Try to open and allocate a buffer
  if( ( tiff = TIFFOpen( filename.c_str(), "rm" ) ) == NULL ){
    throw file_error( "tiff open failed for: " + filename );
  }

  // Load our metadata if not already loaded
  if( bpc == 0 ) loadImageInfo( currentX, currentY );

  // Insist on a tiled image
  if( (tile_width == 0) && (tile_height == 0) ){
    throw file_error( "TIFF image is not tiled" );
  }

  isSet = true;

}


void QPTIFFImage::loadImageInfo( int seq, int ang )
{
  tdir_t current_dir;
  int count;
  uint16 colour, samplesperpixel, bitspersample, sampleformat, compression;
  double sminvaluearr[4] = {0.0}, smaxvaluearr[4] = {0.0};
  double *sminvalue = NULL, *smaxvalue = NULL;
  unsigned int w, h;
  string filename;
  char *tmp = NULL;

  currentX = seq;
  currentY = ang;

  // Get the tile and image sizes
  // Ref:  http://www.libtiff.org/man/TIFFSetField.3t.html
  TIFFGetField( tiff, TIFFTAG_TILEWIDTH, &tile_width );
  TIFFGetField( tiff, TIFFTAG_TILELENGTH, &tile_height );
  TIFFGetField( tiff, TIFFTAG_IMAGEWIDTH, &w );
  TIFFGetField( tiff, TIFFTAG_IMAGELENGTH, &h );
  TIFFGetField( tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel );
  TIFFGetField( tiff, TIFFTAG_BITSPERSAMPLE, &bitspersample );
  TIFFGetField( tiff, TIFFTAG_PHOTOMETRIC, &colour );
  TIFFGetField( tiff, TIFFTAG_SAMPLEFORMAT, &sampleformat );
  TIFFGetField( tiff, TIFFTAG_COMPRESSION, &compression );

  if (w < 2000 && h < 2000){
    tile_width = w;
    tile_height = h;
  }

  if ( samplesperpixel == 1 ){
    // Fluorescence QPTIFF and FUSED_TIFF have N channel directories per resolution.
    // For QPTIFF the N+1st directory is the thumbnail.
    channels = 0;

    unsigned int w1 = w, h1 = h;
    while ( w1 == w && h1 == h ){
      TIFFReadDirectory( tiff );
      TIFFGetField( tiff, TIFFTAG_IMAGEWIDTH, &w1 );
      TIFFGetField( tiff, TIFFTAG_IMAGELENGTH, &h1 );
      channels++;
      if (loglevel >= 2)
        logfile << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << "():: "
                << "Dir w, h " << w1 << ", " << h1 << "; Chan " << channels << endl;
    }
  }
  else {
    // Brightfield QPTIFF have one 3-channel RGB directory per resolution.
    TIFFReadDirectory( tiff );  // For QPTIFF, read the thumbnail.
    channels = 3;
  }

  if (loglevel >= 2)
    logfile << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << "():: "
			<< "TIFFCurrentDirectory " << TIFFCurrentDirectory( tiff )
            << ", Channels: " << channels
            << ", sampleformat = " << sampleformat
            << ", compression = " << compression
            << endl;

  bpc = (unsigned int) bitspersample;
  sampleType = (sampleformat == SAMPLEFORMAT_IEEEFP || compression == COMPRESSION_ADOBE_DEFLATE)
    ? FLOATINGPOINT : FIXEDPOINT;

  // // Check for the no. of resolutions in the pyramidal image
  // current_dir = TIFFCurrentDirectory( tiff );
  // TIFFSetDirectory( tiff, 0 );

  // Store the list of image dimensions available
  image_widths.push_back( w );
  image_heights.push_back( h );

  if (format != QPTIFF) {
    // For non-QPTIFF (e.g., FUSED_TIFF), there was no thumbnail directory, so backup one.
    tdir_t backup_dir;
    backup_dir = TIFFCurrentDirectory( tiff );
    if (backup_dir > 0) --backup_dir;
    TIFFSetDirectory( tiff, backup_dir );
  }

  for( count = 0; TIFFReadDirectory( tiff ); count++ ){
	if (loglevel >= 9)
	  logfile << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << "():: "
			  << "TIFFCurrentDirectory " << TIFFCurrentDirectory( tiff )
			  << endl;

    TIFFGetField( tiff, TIFFTAG_IMAGEWIDTH, &w );
    TIFFGetField( tiff, TIFFTAG_IMAGELENGTH, &h );
    logfile << "QPTIFFImage:: " << "directory " << w << " " << h << endl;
    // ignore downsamples smaller than 2K x 2K
    if (w < 2000 && h < 2000) {
      break;
    }
    image_widths.push_back( w );
    image_heights.push_back( h );

	// Move to the last directory of the current channel.
    current_dir = TIFFCurrentDirectory( tiff );
    TIFFSetDirectory( tiff, current_dir + (channels - 1) );
    if (loglevel >= 2)
      logfile << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << "():: "
              << "Channels: " << count << " w, h " << w << ", " << h << endl;
  }
  // Reset the TIFF directory
  TIFFSetDirectory( tiff, current_dir );

  numResolutions = count+1;

  // Handle various colour spaces
  if( colour == PHOTOMETRIC_CIELAB ) {
    colourspace = CIELAB;
  }
  else if( colour == PHOTOMETRIC_MINISBLACK ){
    colourspace = (bpc==1)? BINARY : GREYSCALE;
  }
  else if( colour == PHOTOMETRIC_PALETTE ){
    // Watch out for colourmapped images. These are stored as 1 sample per pixel,
    // but are decoded to 3 channels by libtiff, so declare them as sRGB
    colourspace = sRGB;
    channels = 3;
  }
  else if( colour == PHOTOMETRIC_YCBCR ){
    // JPEG encoded tiles can be subsampled YCbCr encoded. Ask to decode these to RGB
    TIFFSetField( tiff, TIFFTAG_JPEGCOLORMODE, JPEGCOLORMODE_RGB );
    colourspace = sRGB;
  }
  else if( colour == PHOTOMETRIC_RGB ){
    // RGB stored as a single channel.
    colourspace = sRGB;
  }
  else {
    colourspace = sRGB;
  }

  // Get the max and min values for our data type - required for floats
  // This are usually single values per image, but can also be per channel
  // in libtiff > 4.0.2 via http://www.asmail.be/msg0055458208.html

#ifdef TIFFTAG_PERSAMPLE
  if( channels > 1 ){
    TIFFSetField(tiff, TIFFTAG_PERSAMPLE, PERSAMPLE_MULTI);
    TIFFGetFieldDefaulted( tiff, TIFFTAG_SMINSAMPLEVALUE, &sminvalue );
    TIFFGetFieldDefaulted( tiff, TIFFTAG_SMAXSAMPLEVALUE, &smaxvalue );
    TIFFSetField(tiff, TIFFTAG_PERSAMPLE, PERSAMPLE_MERGED);
    if (!sminvalue) sminvalue = sminvaluearr;
    if (!smaxvalue) smaxvalue = smaxvaluearr;
  }
  else{
#endif
    sminvalue = sminvaluearr;
    smaxvalue = smaxvaluearr;
    TIFFGetFieldDefaulted( tiff, TIFFTAG_SMINSAMPLEVALUE, sminvalue );
    TIFFGetFieldDefaulted( tiff, TIFFTAG_SMAXSAMPLEVALUE, smaxvalue );
#ifdef TIFFTAG_PERSAMPLE
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
  if( TIFFGetField( tiff, TIFFTAG_ARTIST, &tmp ) ) metadata["author"] = tmp;
  if( TIFFGetField( tiff, TIFFTAG_COPYRIGHT, &tmp ) ) metadata["copyright"] = tmp;
  if( TIFFGetField( tiff, TIFFTAG_DATETIME, &tmp ) ) metadata["create-dtm"] = tmp;
  if( TIFFGetField( tiff, TIFFTAG_IMAGEDESCRIPTION, &tmp ) ) metadata["subject"] = tmp;
  if( TIFFGetField( tiff, TIFFTAG_SOFTWARE, &tmp ) ) metadata["app-name"] = tmp;
  if( TIFFGetField( tiff, TIFFTAG_XMLPACKET, &count, &tmp ) ) metadata["xmp"] = string(tmp,count);
  if( TIFFGetField( tiff, TIFFTAG_ICCPROFILE, &count, &tmp ) ) metadata["icc"] = string(tmp,count);

}


void QPTIFFImage::closeImage()
{
  if( tiff != NULL ){
    TIFFClose( tiff );
    tiff = NULL;
  }
  if( tile_buf != NULL ){
    _TIFFfree( tile_buf );
    tile_buf = NULL;
  }
}


RawTile QPTIFFImage::getTile( int seq, int ang, unsigned int res, int layers, unsigned int tile )
{
  uint32 im_width, im_height, tw, th, ntlx, ntly;
  uint32 rem_x, rem_y;
  uint16 colour;
  string filename;


  if (loglevel >= 2)
    logfile << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << "():: "
            << "getTile begin" << endl;

  // Check the resolution exists
  if( res > numResolutions ){
    ostringstream error;
    error << "QPTIFFImage :: Asked for non-existent resolution: " << res;
    throw file_error( error.str() );
  }


  // If we are currently working on a different sequence number, then
  //  close and reload the image.
  if( (currentX != seq) || (currentY != ang) ){
    closeImage();
  }


  // Open the TIFF if it's not already open
  if( !tiff ){
    filename = getFileName( seq, ang );
    if( ( tiff = TIFFOpen( filename.c_str(), "rm" ) ) == NULL ){
      throw file_error( "tiff open failed for:" + filename );
    }
  }


  // Reload our image information in case the tile size etc is different
  if( (currentX != seq) || (currentY != ang) ){
    loadImageInfo( seq, ang );
  }


  // The first resolution is the highest, so we need to invert
  //  the resolution - can avoid this if we store our images with
  //  the smallest image first.
  int vipsres = ( numResolutions - 1 ) - res;

  // Brightfield QPTIFF have one 3-channel RGB directory per resolution.
  // Fluorescence QPTIFF and FUSED_TIFF have N channel directories per resolution.
  // For QPTIFF the N+1st directory is the thumbnail.
  int tiff_dir = vipsres * channels;
  if ( tiff_dir > 0  && format == QPTIFF) {
    tiff_dir += 1;
  }

  if (loglevel >= 2)
    logfile << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << "():: "
                << "tiff_dir = " << tiff_dir << ", bpc = " << bpc << endl;

  // Change to the right directory for the resolution
  if( !TIFFSetDirectory( tiff, tiff_dir ) ) {
    throw file_error( "TIFFSetDirectory failed" );
  }


  // Check that a valid tile number was given
  if( tile >= TIFFNumberOfTiles( tiff ) ) {
    ostringstream tile_no;
    tile_no << "Asked for non-existent tile: " << tile;
    throw file_error( tile_no.str() );
  }


  // Get the size of this tile, the current image,
  //  the number of samples and the colourspace.
  // TIFFTAG_TILEWIDTH give us the values for the resolution,
  //  not for the tile itself
  TIFFGetField( tiff, TIFFTAG_TILEWIDTH, &tw );
  TIFFGetField( tiff, TIFFTAG_TILELENGTH, &th );
  TIFFGetField( tiff, TIFFTAG_IMAGEWIDTH, &im_width );
  TIFFGetField( tiff, TIFFTAG_IMAGELENGTH, &im_height );
  TIFFGetField( tiff, TIFFTAG_PHOTOMETRIC, &colour );
//   TIFFGetField( tiff, TIFFTAG_SAMPLESPERPIXEL, &channels );
//   TIFFGetField( tiff, TIFFTAG_BITSPERSAMPLE, &bpc );

  bool stripped_image = false;

  if (im_width < 2000 && im_height < 2000){
    stripped_image = true;
    tw = im_width;
    th = im_height;
  }
  // Make sure this resolution is tiled
  else if( (tw == 0) || (th == 0) ){
    throw file_error( "Requested resolution is not tiled" );
  }

  // Total number of pixels in tile
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
  if( colour == PHOTOMETRIC_CIELAB ) {
    colourspace = CIELAB;
  }
  else if( colour == PHOTOMETRIC_MINISBLACK ){
    colourspace = (bpc==1)? BINARY : GREYSCALE;
  }
  else if( colour == PHOTOMETRIC_PALETTE ){
    // Watch out for colourmapped images. These are stored as 1 sample per pixel,
    // but are decoded to 3 channels by libtiff, so declare them as sRGB
    colourspace = sRGB;
    channels = 3;
  }
  else if( colour == PHOTOMETRIC_YCBCR ){
    // JPEG encoded tiles can be subsampled YCbCr encoded. Ask to decode these to RGB
    TIFFSetField( tiff, TIFFTAG_JPEGCOLORMODE, JPEGCOLORMODE_RGB );
    colourspace = sRGB;
  }
  else if( colour == PHOTOMETRIC_RGB ){
    // RGB stored as a single channel.
    colourspace = sRGB;
  }
  else {
    colourspace = sRGB;
  }

  tsize_t tile_size = 0;
  if (stripped_image){
    tile_size = im_height * im_height * bpc;
  }
  else {
    tile_size = TIFFTileSize(tiff);
  }

  if (loglevel >= 2)
    logfile << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << "():: "
            << "TIFFTileSize = " << tile_size << endl;

  tdata_t *channels_tile_bufs = new tdata_t[channels];
  for ( int i = 0; i < channels; i++ ) {
    if (loglevel >= 9)
      logfile << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << "():: "
              << "i = " << i << " (of " << channels << ")" << endl;

    if( ( channels_tile_bufs[i] = _TIFFmalloc ( tile_size ) ) == NULL ) {
      throw file_error( "tiff malloc tile failed" );
    }

    if (stripped_image){
      tdata_t buffer = channels_tile_bufs[i];
      int length = 0;
      for (tstrip_t strip = 0; strip < TIFFNumberOfStrips( tiff ); strip++ ){
        length = TIFFReadEncodedStrip( tiff, strip, buffer, (tsize_t) -1 );
        if ( length == -1 ){
          throw file_error( "TIFFReadEncodedStrip failed for " + getFileName( seq, ang ) );
        }
        buffer += length;
      }
    }
    else {
      int length = TIFFReadEncodedTile( tiff, (ttile_t) tile,
              channels_tile_bufs[i], (tsize_t) - 1 );
      if( length == -1 ) {
        if (loglevel >= 2)
          logfile << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << "():: "
                  << "TIFFReadEncodedTile failed for channel i = " << i
                  << " (of " << channels << ")" << endl;
        throw file_error( "TIFFReadEncodedTile failed for " + getFileName( seq, ang ) );
      }

    }
    TIFFReadDirectory( tiff );
  }   

  // Allocate memory for our tile.
  if( !tile_buf ){
    if (loglevel >= 2)
      logfile << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << "():: "
              << "tile_buf size = " << tile_size * channels << endl;

    if( ( tile_buf = _TIFFmalloc( tile_size * channels ) ) == NULL ){
      throw file_error( "tiff malloc tile failed" );
    }
  }

  if (loglevel >= 2)
    logfile << __FILE__ << ":" << __LINE__ << " " << __FUNCTION__ << "():: "
            << "read tiles into channels_tile_bufs"
            << ", bpc = " << bpc << endl;

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
    _TIFFfree(channels_tile_bufs[i]);
  }
  delete channels_tile_bufs;

  RawTile rawtile( tile, res, seq, ang, tw, th, channels, bpc );
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
      // Count backwards as TIFF is usually MSB2LSB
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

}

