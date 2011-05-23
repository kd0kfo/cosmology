
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 28 "physcalc.yacc.ypp"

#include <math.h>  /* For math functions, cos(), sin(), etc.  */
#include <stdio.h>
#include <string.h>
#include <complex>
#include <dirent.h>
#include "symrec.h"
#include "functions.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

  static const char* PROMPT_STRING = ">";
  class DavidException;
  int yylex (void);
  extern void yyrestart(FILE *);
  extern char* yytext;
  void yyerror (char const *);
  void do_ls(const char *path);	
  void print_help(const char *keyword);
  void do_funct(symrec* fnct,  struct calcval *result);
  void do_funct(symrec* fnct, symrec* param, struct calcval *result);
  void do_funct(symrec* fnct, struct calcval param, struct calcval *result);
  void do_funct2(symrec* fnct, symrec* param, symrec* param2, struct calcval *result);
  int update_vars(symrec* var);
  void handle_plane(symrec *rec,double& i, double& j, struct calcval*);
  struct calcval ans;
  void print_complex(struct calcval);
  void print_copyright();

 

/* Line 189 of yacc.c  */
#line 107 "physcalc.yacc.cpp"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NUM = 258,
     VAR = 259,
     FNCT = 260,
     IMAG = 261,
     EXIT = 262,
     COPYRIGHT = 263,
     HELP = 264,
     ANS = 265,
     NEG = 266
   };
#endif
/* Tokens.  */
#define NUM 258
#define VAR 259
#define FNCT 260
#define IMAG 261
#define EXIT 262
#define COPYRIGHT 263
#define HELP 264
#define ANS 265
#define NEG 266




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 60 "physcalc.yacc.ypp"

  struct calcval val;   /* For returning numbers.  */
  struct symrec  *tptr;   /* For returning symbol-table pointers.  */



/* Line 214 of yacc.c  */
#line 172 "physcalc.yacc.cpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 184 "physcalc.yacc.cpp"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   129

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  25
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  4
/* YYNRULES -- Number of rules.  */
#define YYNRULES  28
/* YYNRULES -- Number of states.  */
#define YYNSTATES  57

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   266

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,    22,     2,     2,
      19,    20,    14,    13,    21,    12,     2,    15,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    18,
       2,    11,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    23,     2,    24,    17,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    16
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     6,     7,     9,    12,    18,    21,    24,
      27,    30,    32,    35,    37,    39,    43,    47,    52,    57,
      64,    68,    72,    76,    80,    84,    87,    91,    95
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      26,     0,    -1,    26,    27,    -1,    -1,    18,    -1,     8,
      18,    -1,     9,    19,     5,    20,    18,    -1,     9,    18,
      -1,     7,    18,    -1,    28,    18,    -1,     1,    18,    -1,
       3,    -1,     3,     6,    -1,     4,    -1,    10,    -1,     4,
      11,    28,    -1,     5,    19,    20,    -1,     5,    19,     4,
      20,    -1,     5,    19,    28,    20,    -1,     5,    19,     4,
      21,     4,    20,    -1,    28,    13,    28,    -1,    28,    12,
      28,    -1,    28,    14,    28,    -1,    28,    15,    28,    -1,
      28,    22,    28,    -1,    12,    28,    -1,    28,    17,    28,
      -1,    19,    28,    20,    -1,     4,    23,    28,    21,    28,
      24,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,    76,    76,    77,    81,    82,    83,    84,    85,    86,
      87,    90,    91,    92,    93,    94,    95,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NUM", "VAR", "FNCT", "IMAG", "EXIT",
  "COPYRIGHT", "HELP", "ANS", "'='", "'-'", "'+'", "'*'", "'/'", "NEG",
  "'^'", "';'", "'('", "')'", "','", "'%'", "'['", "']'", "$accept",
  "input", "stmt", "exp", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,    61,    45,    43,    42,    47,   266,    94,    59,    40,
      41,    44,    37,    91,    93
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    25,    26,    26,    27,    27,    27,    27,    27,    27,
      27,    28,    28,    28,    28,    28,    28,    28,    28,    28,
      28,    28,    28,    28,    28,    28,    28,    28,    28
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     0,     1,     2,     5,     2,     2,     2,
       2,     1,     2,     1,     1,     3,     3,     4,     4,     6,
       3,     3,     3,     3,     3,     2,     3,     3,     6
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       3,     0,     1,     0,    11,    13,     0,     0,     0,     0,
      14,     0,     4,     0,     2,     0,    10,    12,     0,     0,
       0,     8,     5,     7,     0,    25,     0,     0,     0,     0,
       0,     0,     9,     0,    15,     0,    13,    16,     0,     0,
      27,    21,    20,    22,    23,    26,    24,     0,    17,     0,
      18,     0,     0,     0,     6,    28,    19
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,    14,    15
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -15
static const yytype_int8 yypact[] =
{
     -15,    23,   -15,   -12,     6,    -8,    -5,     3,     7,   -14,
     -15,    45,   -15,    45,   -15,    59,   -15,   -15,    45,    45,
      34,   -15,   -15,   -15,    35,    12,    70,    45,    45,    45,
      45,    45,   -15,    45,   103,    81,   -10,   -15,    92,    25,
     -15,   107,   107,    12,    12,    12,   103,    45,   -15,    39,
     -15,    29,    46,    31,   -15,   -15,   -15
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -15,   -15,   -15,   -11
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      25,    18,    26,    18,    23,    24,    16,    34,    35,    38,
      48,    49,    17,    19,    20,    19,    41,    42,    43,    44,
      45,    21,    46,     2,     3,    22,     4,     5,     6,    31,
       7,     8,     9,    10,    33,    11,    52,     4,    36,     6,
      39,    12,    13,    53,    10,    51,    11,    54,     4,     5,
       6,    56,     0,    13,    37,    10,     0,    11,    27,    28,
      29,    30,     0,    31,    13,     0,     0,     0,    33,     0,
      55,    27,    28,    29,    30,     0,    31,    32,     0,     0,
       0,    33,    27,    28,    29,    30,     0,    31,     0,     0,
      40,     0,    33,    27,    28,    29,    30,     0,    31,     0,
       0,     0,    47,    33,    27,    28,    29,    30,     0,    31,
       0,     0,    50,     0,    33,    27,    28,    29,    30,     0,
      31,    29,    30,     0,    31,    33,     0,     0,     0,    33
};

static const yytype_int8 yycheck[] =
{
      11,    11,    13,    11,    18,    19,    18,    18,    19,    20,
      20,    21,     6,    23,    19,    23,    27,    28,    29,    30,
      31,    18,    33,     0,     1,    18,     3,     4,     5,    17,
       7,     8,     9,    10,    22,    12,    47,     3,     4,     5,
       5,    18,    19,     4,    10,    20,    12,    18,     3,     4,
       5,    20,    -1,    19,    20,    10,    -1,    12,    12,    13,
      14,    15,    -1,    17,    19,    -1,    -1,    -1,    22,    -1,
      24,    12,    13,    14,    15,    -1,    17,    18,    -1,    -1,
      -1,    22,    12,    13,    14,    15,    -1,    17,    -1,    -1,
      20,    -1,    22,    12,    13,    14,    15,    -1,    17,    -1,
      -1,    -1,    21,    22,    12,    13,    14,    15,    -1,    17,
      -1,    -1,    20,    -1,    22,    12,    13,    14,    15,    -1,
      17,    14,    15,    -1,    17,    22,    -1,    -1,    -1,    22
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    26,     0,     1,     3,     4,     5,     7,     8,     9,
      10,    12,    18,    19,    27,    28,    18,     6,    11,    23,
      19,    18,    18,    18,    19,    28,    28,    12,    13,    14,
      15,    17,    18,    22,    28,    28,     4,    20,    28,     5,
      20,    28,    28,    28,    28,    28,    28,    21,    20,    21,
      20,    20,    28,     4,    18,    24,    20
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{


    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 5:

/* Line 1455 of yacc.c  */
#line 82 "physcalc.yacc.ypp"
    {print_copyright();printf(">");}
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 83 "physcalc.yacc.ypp"
    {print_help((yyvsp[(3) - (5)].tptr)->name);printf(">");}
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 84 "physcalc.yacc.ypp"
    {print_help(NULL);printf(">");}
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 85 "physcalc.yacc.ypp"
    {YYACCEPT;}
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 86 "physcalc.yacc.ypp"
    { print_complex((yyvsp[(1) - (2)].val)); ans = (yyvsp[(1) - (2)].val);}
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 87 "physcalc.yacc.ypp"
    { yyerrok;  printf("%s",PROMPT_STRING);         }
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 90 "physcalc.yacc.ypp"
    { (yyval.val) = (yyvsp[(1) - (1)].val);yylval.val = (yyvsp[(1) - (1)].val);}
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 91 "physcalc.yacc.ypp"
    { yylval.val.re = 0;yylval.val.im = (yyvsp[(1) - (2)].val).re;(yyval.val) = yylval.val;}
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 92 "physcalc.yacc.ypp"
    { (yyval.val).re = (yyvsp[(1) - (1)].tptr)->value.var[0];(yyval.val).im = (yyvsp[(1) - (1)].tptr)->value.var[1];}
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 93 "physcalc.yacc.ypp"
    { (yyval.val) = ans;             }
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 94 "physcalc.yacc.ypp"
    { (yyval.val) = (yyvsp[(3) - (3)].val); (yyvsp[(1) - (3)].tptr)->value.var[0] = (yyvsp[(3) - (3)].val).re;(yyvsp[(1) - (3)].tptr)->value.var[1]= (yyvsp[(3) - (3)].val).im; }
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 95 "physcalc.yacc.ypp"
    {do_funct((yyvsp[(1) - (3)].tptr),&ans); (yyval.val) = ans;}
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 96 "physcalc.yacc.ypp"
    {  do_funct((yyvsp[(1) - (4)].tptr),(yyvsp[(3) - (4)].tptr),&ans);(yyval.val) = ans;}
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 97 "physcalc.yacc.ypp"
    { do_funct((yyvsp[(1) - (4)].tptr),(yyvsp[(3) - (4)].val),&ans); (yyval.val) = ans;}
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 98 "physcalc.yacc.ypp"
    { do_funct2((yyvsp[(1) - (6)].tptr), (yyvsp[(3) - (6)].tptr), (yyvsp[(5) - (6)].tptr),&ans); (yyval.val) = ans;}
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 99 "physcalc.yacc.ypp"
    { (yyval.val).re = (yyvsp[(1) - (3)].val).re + (yyvsp[(3) - (3)].val).re;(yyval.val).im = (yyvsp[(1) - (3)].val).im+(yyvsp[(3) - (3)].val).im;}
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 100 "physcalc.yacc.ypp"
    { (yyval.val).re = (yyvsp[(1) - (3)].val).re - (yyvsp[(3) - (3)].val).re;(yyval.val).im = (yyvsp[(1) - (3)].val).im-(yyvsp[(3) - (3)].val).im;}
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 101 "physcalc.yacc.ypp"
    { (yyval.val).re = (yyvsp[(1) - (3)].val).re * (yyvsp[(3) - (3)].val).re-(yyvsp[(1) - (3)].val).im*(yyvsp[(3) - (3)].val).im;(yyval.val).im = (yyvsp[(1) - (3)].val).re*(yyvsp[(3) - (3)].val).im+(yyvsp[(1) - (3)].val).im*(yyvsp[(3) - (3)].val).re;}
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 102 "physcalc.yacc.ypp"
    { (yyval.val).re = (yyvsp[(1) - (3)].val).re / (yyvsp[(3) - (3)].val).re; (yyval.val).im = (yyvsp[(1) - (3)].val).im/(yyvsp[(3) - (3)].val).re;}
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 103 "physcalc.yacc.ypp"
    { (yyval.val).re = (int)(yyvsp[(1) - (3)].val).re % (int)(yyvsp[(3) - (3)].val).re;}
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 104 "physcalc.yacc.ypp"
    { (yyval.val).re = -(yyvsp[(2) - (2)].val).re;(yyval.val).im = -(yyvsp[(2) - (2)].val).im;         }
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 105 "physcalc.yacc.ypp"
    { (yyval.val).re = pow ((yyvsp[(1) - (3)].val).re, (yyvsp[(3) - (3)].val).re);(yyval.val).im=0; }
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 106 "physcalc.yacc.ypp"
    { (yyval.val) = (yyvsp[(2) - (3)].val);                         }
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 107 "physcalc.yacc.ypp"
    {handle_plane((yyvsp[(1) - (6)].tptr),(yyvsp[(3) - (6)].val).re,(yyvsp[(5) - (6)].val).re,&ans);(yyval.val) = ans;}
    break;



/* Line 1455 of yacc.c  */
#line 1591 "physcalc.yacc.cpp"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 1675 of yacc.c  */
#line 110 "physcalc.yacc.ypp"


symrec *
     putsym (char const *sym_name, int sym_type)
     {
       symrec *ptr;
       ptr = (symrec *) malloc (sizeof (symrec));
       ptr->name = (char *) malloc (strlen (sym_name) + 1);
       strcpy (ptr->name,sym_name);
       ptr->type = sym_type;
       ptr->value.var[0] = 0; /* Set value to 0 even if fctn.  */
       ptr->value.var[1] = 0;
       ptr->next = (struct symrec *)sym_table;
       sym_table = ptr;
       return ptr;
     }
     
     symrec *
     getsym (char const *sym_name)
     {
       symrec *ptr;
       for (ptr = sym_table; ptr != (symrec *) 0;
            ptr = (symrec *)ptr->next)
         if (strcmp (ptr->name,sym_name) == 0)
           return ptr;
       return 0;
     }

#if 0
#include <ctype.h>

     int
     yylex (void)
     {
       int c;

       /* Ignore white space, get first nonwhite character.  */
       while ((c = getchar ()) == ' ' || c == '\t');
     
       if (c == EOF)
         return 0;
     
       /* Char starts a number => parse the number.         */
       if (c == '.' || isdigit (c))
         {
           ungetc (c, stdin);
           scanf ("%lf", &yylval.val[0]);
           return NUM;
         }
     
       /* Char starts an identifier => read the name.       */
       if (isalpha (c))
         {
           symrec *s;
           static char *symbuf = 0;
           static int length = 0;
           int i;
     
           /* Initially make the buffer long enough
              for a 40-character symbol name.  */
           if (length == 0)
             length = 40, symbuf = (char *)malloc (length + 1);
     
           i = 0;
           do
             {
               /* If buffer is full, make it bigger.        */
               if (i == length)
                 {
                   length *= 2;
                   symbuf = (char *) realloc (symbuf, length + 1);
                 }
               /* Add this character to the buffer.         */
               symbuf[i++] = c;
               /* Get another character.                    */
               c = getchar ();
             }
           while (isalnum (c));
     
           ungetc (c, stdin);
           symbuf[i] = '\0';
     
           s = getsym (symbuf);
           if (s == 0)
             s = putsym (symbuf, VAR);
           yylval.tptr = s;
           return s->type;
         }
     
       /* Any other character is a token by itself.        */
       return c;
     }
#endif

/* Called by yyparse on error.  */
     void
     yyerror (char const *s)
     {
       printf ("%s: %s\n", s,yytext);
     }
     
     struct init
     {
       char const *fname;
       double (*fnct) (double);
       void (*v_c_fnct)(const char*);
     };
     
     struct init const arith_fncts[] =
     {
       "abs",  fabs,0,
       "acos", acos,0,
       "asin", asin,0,
       "atan", atan,0,
       "ceil", ceil,0,
       "cos",  cos,0,
       "cosh",  cosh,0,
       "exp",  exp,0,
       "floor", floor,0,
       "sin",  sin,0,
       "sinh",  sinh,0,
       "sqrt", sqrt,0,
       "tan",  tan,0,
       "tanh",  tanh,0,
       "log",  log10,0,
       "ln",   log,0,
       0, 0, 0
     };

	struct init const other_fncts[] = 
	{
		"ls",0,do_ls,
		0,0,0
	};
     
     /* The symbol table: a chain of `struct symrec'.  */
     symrec *sym_table;
     
/**
 * Initialize function table
 */
     void
     init_table (void)
     {
       int i;
       symrec *ptr;
       const struct init *curr = arith_fncts;

	i = 0;
       while (curr[i].fname != 0)
         {
           ptr = putsym (curr[i].fname, FNCT);
	   if(curr[i].fnct != 0)
	           ptr->value.fnctptr =	 curr[i].fnct;
	   else if(curr[i].v_c_fnct != 0)
	          ptr->value.v_c_fnctptr =	 curr[i].v_c_fnct;
           ptr->plane_fnctptr = NULL;
	

           i++;
	   if(curr == arith_fncts && curr[i].fname == 0)
	   {
		curr = other_fncts;i = 0;
	   }
         }

       putsym("print",FNCT)->plane_fnctptr = print_plane;
       putsym("clear",FNCT)->plane_fnctptr = clear_plane;
       putsym("add",FNCT)->plane_fnctptr = add_planes;
       putsym("subtract",FNCT)->plane_fnctptr = subtract_planes;
       putsym("multiply",FNCT)->plane_fnctptr = multiply_planes;
       putsym("open",FNCT)->plane_fnctptr = open_plane;
       putsym("save",FNCT)->plane_fnctptr = save_plane;
       putsym("copy",FNCT)->plane_fnctptr = copy_plane;
       putsym("fourier",FNCT)->plane_fnctptr = fourier_plane;
       putsym("ifourier",FNCT)->plane_fnctptr = ifourier_plane;
     }
     
     int
     main (int argc, char **argv)
     {
       symrec *funct;
       FILE *input = NULL;
       print_copyright();
       init_table ();
       if(argc > 1)
       {
	 input = fopen(argv[1],"r");
	 if(input != NULL)
	 	  yyrestart(input);
       }

       if(input == NULL)
       {
              printf("\n");
       	      printf("Welcome! To exit, type quit.\nALL statements must end with a semicolon.\n%s",PROMPT_STRING);
	      yyparse();
	      printf("Good bye.\n");
      }
      else
      {
	yyparse();      
      }
     return 0;
}

void do_funct(symrec *rec, struct calcval param, struct calcval *result)
{
  if(rec == NULL || rec->plane_fnctptr != NULL)
    return;

    result->re = (*rec->value.fnctptr)(param.re);
    result->im = 0;
}

void do_funct(symrec *rec, symrec *param, struct calcval *result)
{
  if(rec == NULL || param == NULL)
    return;

  if(rec->plane_fnctptr != NULL)
    {
      symrec* newrec = (*rec->plane_fnctptr)(&param,1);
      if(newrec != NULL)
	{
	  newrec->next = (struct symrec *)sym_table;
	  sym_table = newrec;
	}
       return;
    }

  if(rec->value.v_c_fnctptr != 0)
  {
	(*rec->value.v_c_fnctptr)(param->name);
	return;
  }

  result->re = (*(rec->value.fnctptr))(param->value.var[0]);
  result->im = 0.0;
}


void do_funct2(symrec* rec, symrec* param, symrec* param2, struct calcval *result)
{
  if(rec == NULL || param == NULL || param2 == NULL)
    return;
  if(rec->plane_fnctptr == NULL)
    return;

  symrec* params[] = {param,param2};
  symrec* newrec = (*rec->plane_fnctptr)(params,2);
  if(newrec != NULL)
    {
      newrec->next = (struct symrec *)sym_table;
      sym_table = newrec;
    }

}

void handle_plane(symrec *rec,double& i, double& j, struct calcval *result)
{
  rec->isPlane = true;
  if(rec->value.planeptr == NULL)
    {
      rec->value.planeptr = create_plane(i,j); 
      return;
    }
  try{
    const Double& complex = rec->value.planeptr->getValue((int)i,(int)j);
    result->re = complex.getValue(0);
    result->im = complex.getValue(1);
  }
  catch(DavidException de)
    {
      std::cout << de.what() << std::endl;
    }
}

void print_complex(struct calcval cnumber)
{
  using namespace std;
  cout << "\tans: ";
  if(cnumber.im != 0)
    {
      math::Complex buff(cnumber.re,cnumber.im);
      cout << buff;
    }
  else
    cout << cnumber.re;
  cout <<   endl << ">";
}

void do_funct(symrec* fnct,  struct calcval *result)
{
	if(fnct == NULL)
		return;
	if(fnct->type != FNCT)
		return;
	if(fnct->value.v_c_fnctptr == 0)
		return;
	(*fnct->value.v_c_fnctptr)(NULL);
}

void print_help(const char *keyword)
{
	if(keyword == NULL || strcmp(keyword,"print") == 0)
		printf("print -- summarizes data contined in the specified variable.\n");
	if(keyword == NULL || strcmp(keyword,"clear") == 0)
		printf("clear -- clears plane variable contents. \n");
	if(keyword == NULL || strcmp(keyword,"add") == 0)
		printf("add -- adds two planes. \n\tUsage: add(A,B)\n\t\t A <- A+B\n");
	if(keyword == NULL || strcmp(keyword,"substract") == 0)
		printf("substracts -- subtracts two planes. \n\tUsage: subtract(A,B)\n\t\t A <- A-B\n");
	if(keyword == NULL || strcmp(keyword,"multiply") == 0)
		printf("add -- multiply two planes. \n\tUsage: multiply(A,B)\n\t\t A <- A*B\n");
	if(keyword == NULL || strcmp(keyword,"open") == 0)
		printf("open -- opens a plane and stores it in the specified variable.\n\tUsage: open(A,filename);\n");
	if(keyword == NULL || strcmp(keyword,"save") == 0)
		printf("save -- saves a plane save the specified file.\n\tUsage: save(A,filename);\n");
	if(keyword == NULL || strcmp(keyword,"copy") == 0)
		printf("copy -- copies a plane.\n\tUsage: copy(A,B)\n\t\t B <- A\n");
	if(keyword == NULL || strcmp(keyword,"fourier") == 0)
		printf("fourier -- fourier transforms a plane. Performed in place.\nUsage: fourier(A)\n\t\tA <- F(A)\n");
	if(keyword == NULL || strcmp(keyword,"ifourier") == 0)
		printf("ifourier -- inverse fourier transforms a plane. Performed in place.\nUsage: ifourier(A)\n\t\tA <- F^-1 (A)\n");
	if(keyword == NULL)
	{
		printf("\nThe following standard functions are also included:\n"); 
		const struct init *curr = arith_fncts;
		while(curr != NULL && curr->fname != 0)
		{
			printf("%s",curr->fname);
			curr = curr++;
			if(curr != NULL && curr->fname != 0)
				printf(",");
		}
		printf("\n");
	}
}

void do_ls(const char *path)
{
	DIR *dp;
	struct dirent *ep;
	
	if(path == NULL)
		dp = opendir ("./");
	else
		dp = opendir(path);
	if (dp != NULL)
	{
		while (ep = readdir (dp))
			if(ep->d_name && ep->d_name[0] != '.')
				puts (ep->d_name);
		closedir (dp);
	}
	else
		printf ("Couldn't open the directory\n");

	
}

void print_copyright()
{
	printf("physcalc v%s  Copyright (C) 2010 David Coss, PhD\n",PACKAGE_VERSION);
	printf("This program comes with ABSOLUTELY NO WARRANTY; for details visit http://www.gnu.org/licenses/gpl-2.0.html or read the COPYING file distributed with this program.\n");
}
