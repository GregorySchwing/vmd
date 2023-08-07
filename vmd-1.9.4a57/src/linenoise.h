/*
 * $Id: linenoise.h,v 1.3 2020/05/22 05:53:27 johns Exp $
 *
 * Modified linenoise from the standard distribution by adding 
 * user context parameters for each of the callbacks, so that
 * the caller can pass in runtime-generated lists of command
 * completion strings without the use of global variables or 
 * other undesirable methods.
 *
 * Further modifications were made to the linenoise Unix TTY handling code
 * to allow VMD to control TTY buffer flushes when switching between raw and 
 * cooked TTY mode.  This is required so that VMD can switch to raw mode
 * for character-at-a-time input, so that linenoise doesn't have 
 * blocking behavior that prevents the VMD main loop from free-running. 
 * With these changes, VMD free runs except when in actual command editing.  
 * Correct handling of VMD console output is made more complex by entry 
 * and return from raw TTY mode, but it is a much more usable scenario
 * for the end user.
 */

/* linenoise.h -- VERSION 1.0
 *
 * Guerrilla line editing library against the idea that a line editing lib
 * needs to be 20,000 lines of C code.
 *
 * See linenoise.c for more information.
 *
 * ------------------------------------------------------------------------
 *
 * Copyright (c) 2010-2014, Salvatore Sanfilippo <antirez at gmail dot com>
 * Copyright (c) 2010-2013, Pieter Noordhuis <pcnoordhuis at gmail dot com>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *  *  Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *  *  Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __LINENOISE_H
#define __LINENOISE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct linenoiseCompletions {
  size_t len;
  char **cvec;
} linenoiseCompletions;

typedef void(linenoiseCompletionCallback)(void *, const char *, linenoiseCompletions *);
typedef char*(linenoiseHintsCallback)(void *, const char *, int *color, int *bold);
typedef void(linenoiseFreeHintsCallback)(void *, void *);
void linenoiseSetCompletionCallback(void *, linenoiseCompletionCallback *);
void linenoiseSetHintsCallback(void *, linenoiseHintsCallback *);
void linenoiseSetFreeHintsCallback(void *, linenoiseFreeHintsCallback *);
void linenoiseAddCompletion(linenoiseCompletions *, const char *);

char *linenoise(const char *prompt);
void linenoiseFree(void *ptr);
int linenoiseHistoryAdd(const char *line);
int linenoiseHistorySetMaxLen(int len);
int linenoiseHistorySave(const char *filename);
int linenoiseHistoryLoad(const char *filename);
void linenoiseClearScreen(void);
void linenoiseSetMultiLine(int ml);
void linenoisePrintKeyCodes(void);
void linenoiseMaskModeEnable(void);
void linenoiseMaskModeDisable(void);

int enableRawMode(int fd, int flush);
void disableRawMode(int fd, int flush);

#ifdef __cplusplus
}
#endif

#endif /* __LINENOISE_H */
