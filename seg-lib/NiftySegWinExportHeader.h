#ifndef __NiftySegWinExportHeader_h
#define __NiftySegWinExportHeader_h

#if (defined(_WIN32) || defined(WIN32) || defined(_WINDOWS)) && !defined(LINK_STATIC)
  #ifdef NIFTYSEG_WINDOWS_EXPORT
    #define NIFTYSEG_WINEXPORT __declspec(dllexport)
  #else
    #define NIFTYSEG_WINEXPORT
  #endif  /* NIFTYSEG_WINEXPORT */
#else
/* linux/mac needs nothing */
  #define NIFTYSEG_WINEXPORT
#endif

#endif  //NiftySegWinExportHeader
