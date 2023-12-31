/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2019 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/
/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: win32vmdstart.c,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.52 $      $Date: 2020/12/14 20:26:02 $
 *
 ***************************************************************************/
/**
 *  \file win32vmdstart.c
 *  \brief Win32/64 platform specific application startup code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <windows.h>

#define VMD_FILENAME_MAX 1024

static char * vmd_get_vmddir(void) {
  LONG res = ERROR_FILE_NOT_FOUND;
  HKEY k;
  BYTE * vmddir;
  DWORD bufsz;
  DWORD buftype;

#if defined(_WIN64)
  const char *vmdregkey = "Software\\University of Illinois\\VMD\\1.9.4";
#else
  const char *vmdregkey = "Software\\University of Illinois\\VMD\\x86\\1.9.4";
#endif
  res = RegOpenKeyEx(HKEY_LOCAL_MACHINE, vmdregkey, 0, KEY_READ, &k);
  if (res != ERROR_SUCCESS)  {
    printf("Failed to find registry key: '%s'.\n", vmdregkey);
    return NULL;
  }

  bufsz = 2048;
  vmddir = (BYTE *) malloc(2048);

  res = RegQueryValueEx(k, "VMDDIR", NULL, &buftype, vmddir, &bufsz);

  if (res != ERROR_SUCCESS) {
    printf("Failed to get VMDDIR registry key value.\n");
    free(vmddir);
    vmddir = NULL;
  }

  if (buftype != REG_SZ) {
    printf("Bad data type in VMDDIR registry key value\n");
    free(vmddir);
    vmddir = NULL;
  }

  return (char *) vmddir;
}

#if defined(VMDSEPARATESTARTUP)
// used when seperate startup (for shell) is needed
int main(int argc, char **argv) {
#else
// called by VMD main program
int win32vmdstart(void) {
#endif

  int i, j, len;
  char *str, *str2, *vmddir;
  char tmp[MAX_PATH + 8192];
  char buffer[MAX_PATH +1];
  int usedregistry;
#if defined(VMDSEPARATESTARTUP)
  int rc;
#endif

/* XXX debugging hacks to be sure the compiler mode is set correctly */
#if 0 && defined(_WIN64)
  printf("TEST: Windows x64 VMD test build.\n");
  printf("TEST: sizeof(int): %d  sizeof(long): %d  sizeof(void *): %d\n", 
    sizeof(int), sizeof(long), sizeof(void *));
#endif 

  vmddir = NULL;             // initialize to NULL for safety
  usedregistry = 1;          // change this to 0 if we fail
  vmddir = vmd_get_vmddir(); // search Windows registry first 

  if (vmddir == NULL) {
#if 1
    // get full pathname to VMD executable
    GetModuleFileName(NULL, buffer, sizeof(buffer));
    usedregistry = 0;  // registry failed, we had to guess, using exe path
    vmddir = buffer;   // get from the Win32 API vmd.exe path
#elif defined(VMDSEPARATESTARTUP)
    // desparate attempt to get VMDDIR from argv[0], which only works if 
    // we were started from the explorer shell
    usedregistry = 0;  // registry failed, we had to guess, using exe path
    vmddir = argv[0];  // get from the explorer shell vmd.exe path
#else
    return -1; // fail and exit
#endif
  }

  len = strlen(vmddir);
  str = (char *) malloc(len);
  str2 = (char *) malloc(len);
  strcpy(str,  vmddir);
  strcpy(str2, vmddir);

  j=len;
  for (i=0; i<len; i++) {
    if (str[i] == '\\') {
      str2[i] = '/';
      j=i;
    }
  }

  if (usedregistry) {
    free(vmddir);  // free memory allocated by the vmd_get_vmddir() function
    vmddir = NULL; // make sure we don't accidentally use it again
  } else {  
    // If we didn't use the registry, we need to strip off everything after
    // and including the last "\" character. 
    // If we did use the registry, no string truncation is required.
    str[j] = '\0';
    str2[j] = '\0';
  }

  strcpy(tmp, "VMDDIR=");
  strcat(tmp, str2);
  putenv(strdup(tmp));

  strcpy(tmp, "TCL_LIBRARY=");
  strcat(tmp, str2);
  strcat(tmp, "/scripts/tcl");
  putenv(strdup(tmp));

  strcpy(tmp, "TK_LIBRARY=");
  strcat(tmp, str2);
  strcat(tmp, "/scripts/tk");
  putenv(strdup(tmp));

  if (!getenv("PYTHONPATH")) {
    strcpy(tmp, "PYTHONPATH=");
    strcat(tmp, str2);
    strcat(tmp, "/scripts/python");
    putenv(strdup(tmp));
  } else {
    strcpy(tmp, getenv("PYTHONPATH"));
    strcat(tmp, ":");
    strcat(tmp, str2);
    strcat(tmp, "/scripts/python");
    putenv(strdup(tmp));
  }

  strcpy(tmp, "VMDBABELBIN=");
  strcat(tmp, str);
  strcat(tmp, "\\babel\\babel.exe");
  putenv(strdup(tmp));

  strcpy(tmp, "BABEL_DIR=");
  strcat(tmp, str);
  strcat(tmp, "\\babel");
  putenv(strdup(tmp));

#if defined(_WIN64)
  strcpy(tmp, "STRIDE_BIN=\"");
  strcat(tmp, str);
  strcat(tmp, "\\stride_WIN64.exe\"");
  putenv(strdup(tmp));

  strcpy(tmp, "SURF_BIN=");
  strcat(tmp, str);
  strcat(tmp, "\\surf_WIN64.exe");
  putenv(strdup(tmp));

  strcpy(tmp, "TACHYON_BIN=");
  strcat(tmp, str);
  strcat(tmp, "\\tachyon_WIN64.exe");
  putenv(strdup(tmp));
#else
  strcpy(tmp, "STRIDE_BIN=\"");
  strcat(tmp, str);
  strcat(tmp, "\\stride_WIN32.exe\"");
  putenv(strdup(tmp));

  strcpy(tmp, "SURF_BIN=");
  strcat(tmp, str);
  strcat(tmp, "\\surf_WIN32.exe");
  putenv(strdup(tmp));

  strcpy(tmp, "TACHYON_BIN=");
  strcat(tmp, str);
  strcat(tmp, "\\tachyon_WIN32.exe");
  putenv(strdup(tmp));
#endif

#if defined(VMDSEPARATESTARTUP)
  strcpy(tmp, str);
  strcat(tmp, "\\winvmd.exe");

  // Add all of the command line options we may have had to the 
  // invocation of the real program.
  strcat(tmp, " ");
  for(i=1; i<argc; i++) {
    strcat(tmp, argv[i]);
    strcat(tmp, " ");
  }

  // Run the program!
  rc = WinExec(tmp, SW_SHOW);
  if (!rc) {
    printf("failed to start program %s\n", tmp);
    Sleep(10000);  
  }
#endif

  return 0;
}




