#ifndef PHYSCALC_H
#define PHYSCALC_H 1
#include "symrec.h"

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
