/******************************************************************************
　　　　　　　　simplex.h

目的：simplex_io.c および simplex.c のためのヘッダファイル
						   作成日　2000.8.24-2000.9.18
						   修正日  2001.10.16
******************************************************************************/

#include<stdio.h>               /* standard input output */
#include<stdlib.h>              /* standard library */
#include<string.h>              /* 文字列を操作するためのヘッダファイル */
#include<assert.h>		/* 診断機能を付加するためのヘッダファイル */

#define CONSTRAIN_MAX 50	/* 制約条件式の最大数 */
#define VARIABLE_MAX 90		/* 余裕変数やスラック変数を含む変数の最大数 */
#define MAXWORD 80		/* 読み込みファイル１行の最大文字数 */
#define INF 99999.0             /* 無限大 */
#define ZI 0.001                /* ゼロ判定値 */

enum maxmin_type {max,min};	/* 列挙型の宣言 */
enum yes_no{yes,no};

/*-----------------------------------------------------------------------------
                ファイル読み込み用 構造体
-----------------------------------------------------------------------------*/

typedef struct
{
	int x_sub;		/* Ｘ１、Ｘ２などＸの添字（一時記憶用） */
	double coef;		/* 係数 */	
	char line[MAXWORD];     /* １行分の文字列を格納する配列 */
	char *string;	        /* 文字配列ポインタ */

} input_file;                   /* input_file型の構造体としてつけた名前 */

/*-----------------------------------------------------------------------------
		プログラム条件に関するグローバル変数
-----------------------------------------------------------------------------*/

int m;	             /* 制約条件式の数 */
int n;	             /* 変数Ｘの数 */
int accuracy;	     /* 計算精度 */
int nn;	             /* 変数の総数(Ｘと余裕変数・スラック変数の数の合計) */

/*-----------------------------------------------------------------------------
                目的関数に関するグローバル変数

	z = c[1] * x[1] + c[2] * x[2] + ‥‥‥ + c[n] * x[n]
-----------------------------------------------------------------------------*/

double z;	                        /* 目的関数の値 */
double c[VARIABLE_MAX+CONSTRAIN_MAX];	/* 目的関数の係数 */
enum maxmin_type max_min;		/* maxmin_type型の変数max_minの宣言 */

/*-----------------------------------------------------------------------------
                制約条件に関するグローバル変数

	a[1][1] * x1 + a[1][2] * x2 + ‥‥‥ + a[1][n] * xn sign[1] b[1]
	a[2][1] * x1 + a[2][2] * x2 + ‥‥‥ + a[2][n] * xn sign[2] b[1]
	‥‥‥‥‥‥‥‥‥
	a[m][1] * x1 + a[m][2] * x2 + ‥‥‥ + a[m][n] * xn sign[m] b[1]
-----------------------------------------------------------------------------*/

double a[CONSTRAIN_MAX][VARIABLE_MAX+CONSTRAIN_MAX];
		/* 制約条件式の係数（余裕変数・スラック変数も含む）*/
double b[CONSTRAIN_MAX];		/* 制約条件式の右辺定数 */
char sign[CONSTRAIN_MAX];		/* 等号不等号 */

/*-----------------------------------------------------------------------------
                演算に関するグローバル変数
-----------------------------------------------------------------------------*/

int artivari[VARIABLE_MAX+CONSTRAIN_MAX];    /* 人為変数 */
int artivari_num;                            /* 人為変数の数 */
int basic_variable[CONSTRAIN_MAX];	     /* 基底変数となる変数の番号 */

/*-----------------------------------------------------------------------------
                関数のプロトタイプ宣言 (名前と引数だけの原型宣言)
-----------------------------------------------------------------------------*/
/* ファイル:simplex_io.c */
void read_data();
void read_obj_func(input_file*, FILE *);    /* FILE*はファイルポインタ */
void read_constrain(input_file*, FILE *);
void read_one_line(FILE *,char *);
int read_one_token(char *,double *);
void output_title(FILE *);
/* ファイル:simplex.c */
void calculate_simplex(FILE *);
void init();
void simplex_criterion();
/* 第１フェーズ */
void reset_omega(double *);    /* double型変数を示すポインタ */        
int get_first(double *);
int get_max_omega(double *);
void calculate_theta(double *,int);
int get_min_theta(double *);
void base_change(int,int,double *);
void pivoting(int,int,double,double *);    /* doubleはdouble型の変数 */
/* 第２フェーズ */
int get_second();
int get_max_c();
void reset_theta(double *);
/* 出力 */
void output_tableau(FILE *,double *);
void output_omega(FILE *,double *);
void output_result(FILE *);
double result_read(int);     /* intは整数型の何らかの変数 */

/* doubleは単に値を参照するための引数である */
/* double*は値を参照し、かつ値を変更して返り値とするための引数 */
