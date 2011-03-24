%{
#include <stdlib.h>
#include "mycosmo.tab.h"
#include "symrec.h"

void yyerror(const char *);


%}

%% 
"ans"  return ANS;
"exit" return EXIT;
"quit" return EXIT;

[a-zA-Z]*  {
           symrec* rec = getsym(yytext);
           if(rec == NULL)
               rec = putsym(yytext,VAR);
           yylval.tptr = rec;
          return rec->type;
        }

[0-9]*    {
            yylval.val = atoi(yytext);
            return NUM;
          }

[-()%<>=+*/;\[\],{}.^]  {  return *yytext;  }

[ \t\n]+  ;//ignore whitespace


.  yyerror("Unknown Character.");

%%

int yywrap(void)
{
printf(">");
return 1;
}
