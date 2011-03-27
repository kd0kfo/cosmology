%{
#include <stdlib.h>
#include "physcalc.tab.h"
#include "symrec.h"

#ifndef M_PI
#define M_PI 3.1415926
#endif

void yyerror(const char *);


%}

%% 
"ans"  return ANS;
"exit" return EXIT;
"quit" return EXIT;
"i" return IMAG;

"pi" {yylval.val[0] = M_PI;yylval.val[1] = 0; return NUM;}

[a-zA-Z][a-zA-Z0-9]*  {
           symrec* rec = getsym(yytext);
           if(rec == NULL)
               rec = putsym(yytext,VAR);
           yylval.tptr = rec;
          return rec->type;
        }

[0-9]*\.[0-9]+ {
                  yylval.val[0] = atof(yytext);yylval.val[1] = 0;
                  return NUM;
              }

[0-9]*    {
            yylval.val[0] = atoi(yytext);yylval.val[1] = 0.0;
            return NUM;
          }

[-()%<>=+*/;\[\],{}^]  {  return *yytext;  }

[ \t\n]+  ;//ignore whitespace


.  yyerror("Unknown Character.");

%%

int yywrap(void)
{
printf(">");
return 1;
}

