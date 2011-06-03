/**
 * 
 * This file is part of physcalc, an interactive calculator utility
 * designed to carry out lensing calculations and to manipulate plane
 * data structures.
 *
 * This file contains parser code. It is meant to be used with
 * yacc (specificly GNU bison). For an excellent "A Compact Guide
 * to Lex & Yacc" by Tom Niemann" available at epaperpress.com. 
 * Portions of source code from that paper are used here with 
 * permission given in the Preface (page 3).
 * 
 * Copyright 2007, 2010 David Coss, PhD
 *
 * physcalc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * physcalc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with physcalc.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef PHYSCALC_H
#define PHYSCALC_H 1
#include "symrec.h"

extern struct calcval PHYSCALC_ans;

extern const char* PROMPT_STRING;
extern void yyrestart(FILE *);
extern char* yytext;
void do_ls(const char *path);	
void print_help(const char *keyword);
void do_funct(symrec* fnct,  struct calcval *result);
void do_funct(symrec* fnct, symrec* param, struct calcval *result);
void do_funct(symrec* fnct, struct calcval param, struct calcval *result);
void do_funct2(symrec* fnct, symrec* param, symrec* param2, struct calcval *result);
int update_vars(symrec* var);
void handle_plane(symrec *rec,double& i, double& j, struct calcval*);

void print_complex(struct calcval);
void print_copyright();
void init_table();

#endif
