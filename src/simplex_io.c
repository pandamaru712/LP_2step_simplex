/******************************************************************************
		ファイル名：simplex_io.c：入出力のためのファイル
		内容：ファイルから入力データを読み込み,結果を出力する
******************************************************************************/

#include "simplex.h"

/*-----------------------------------------------------------------------------
　　　関数名：read_data
	　内容：ファイルからデータを読み込む関数
	　引数：なし
-----------------------------------------------------------------------------*/
void read_data()
{
  /***** 初期設定 **********************************************************/
  int i,j;		     /* カウンタ */
  char filename[MAXWORD];    /* 入力データのファイルネーム */
  input_file input;	     /*input_file型の変数inputを定義する*/
  FILE *fp;		     /* ファイルポインタ */
  
  /***** ファイルの読み込み ***********************************************/
  printf("問題式のファイルネームを入力して下さい\nFileName?");
  scanf("%s",filename);
  if((fp=fopen(filename,"r"))==NULL)
    {	/* ファイルがオープンできなかったとき */
      printf("File Open Error!! オープンできません\n");
      exit(0);
    }
  else
    printf("%sファイルのオープンに成功\n",filename);
  /***** プログラム条件の読み込み *****/
  printf("\nプログラム条件読み込み中・・・\n");
  fscanf(fp,"%d %d %d ",&m, &n, &accuracy);  /* &はアドレス演算子 */
  printf("制約条件式 %d個\n",m);
  printf("変数Ｘの数 %d個\n",n);
  printf("計算精度   %d\n",accuracy);
  
  /***** 配列の初期化 *****************/
  for(i=1;i <= m;i++){
    for(j=1;j <= n;j++){
      a[i][j]=0.0;	
    }
    basic_variable[i] = -1*i;
  }
  
  /***** 目的関数の読み込み ***********/
  printf("\n目的関数読み込み中・・・\n");
  read_obj_func(&input, fp);
  
  /***** 制約条件式の読み込み *********/
  printf("\n制約条件式読み込み中・・・\n");
  read_constrain(&input, fp);
  
  fclose(fp);
}
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		関数名：read_obj_func 
		内容：目的関数の読み込み
		引数：input_file & FILE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void read_obj_func(input_file *input, FILE *fp)
{
  int i;
  read_one_line(fp,input->line);    /* ファイルからデータを１行読み込む */
  /* 目的関数の’Ｚ’および’＝’を読み込む */
  input->string = strtok(input->line," ");
  assert(input->string[0]=='z'||input->string[0]=='Z');  /* ->はアロー演算子 */
  /* １文字目がZでなければプログラムを終了する */
  printf("%c ",input->string[0]);
  input->string=strtok(NULL," ");   /* 引き続きスペースの前まで読み込む */
  assert(input->string[0]=='=');	
  /* 次の１文字目が＝でなければプログラムを終了する */
  printf("%c ",input->string[0]);
  /* 目的関数の係数を読み込む */
  for(i=1; i<=n; i++)		
    {
      input->string=strtok(NULL," ");		 /* トークンを読み込む */
      input->x_sub = read_one_token(input->string, &input->coef);
      /* 読み込んだトークンを構文解析する */
      c[input->x_sub]= input->coef;					
      /* 目的関数の係数Cに代入する */
      printf("%6.2fX%d ",input->coef,input->x_sub);
    }
  /* 目的関数の符号（最大化or最小化）を読み込む*/
  input->string=strtok(NULL," ");
  if((input->string[1]=='a')||(input->string[1]=='A'))
    {
      max_min = max;
      printf(" MAX\n");
    }
  else if((input->string[1]=='i')||(input->string[1]=='I'))
    {
      max_min = min;
      printf(" MIN\n");
    }
  else
    { /* エラーコードの出力 */
      printf("目的関数の最大化or最小化のデータが不正です!!\n");
      assert((input->string[1]=='a')&&(input->string[1]=='A')&&
	     (input->string[1]=='i')&&(input->string[1]=='I'));
    }
}
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		関数名：read_constrain
		内容：制約条件式の読み込み
		引数：input_file & FILE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void read_constrain(input_file *input, FILE *fp)
{
  int i;
  for(i=1; i<=m; i++)		        /* 制約条件式の数だけ繰り返す */
    {
      printf("第%2d式  ",i);
      read_one_line(fp,input->line);    /* ファイルからデータを１行読み込む */
/*--------------------------------------
	制約条件式の解析 
--------------------------------------*/
      /* 1:左辺値の読み込み */
      input->string=strtok(input->line," ");	
      /* 最初は先頭からトークンを読み込む */
      while((input->string[0]!='<')&&(input->string[0]!='=')&&
	    (input->string[0]!='>')&&(input->string[0]!='\n'))
	{ /* トークンが不等号か改行コードを含まない間繰り返す */
	  input->x_sub = read_one_token(input->string, &input->coef);
	  /* 読み込んだトークンを構文解析する */
	  a[i][input->x_sub] = input->coef;   
	  /* 制約条件式の係数aに代入する */
	  printf("%6.2fX%d ",input->coef,input->x_sub);
	  input->string=strtok(NULL," ");	/* 新しいトークンを読み込む */
	}
      /* 2:不等号の読み込み */ 
      sign[i]=input->string[0];
      printf("%c= ",sign[i]);
      if(input->string[0]=='\n')
	{ /* エラーチェック（不等号でなく改行コードだったとき） */
	  printf("データに正しい不等号がありません\n");
	  assert(input->string[0]!='\n');
	}
      /* 3:右辺の定数の読み込み */
      input->string=strtok(NULL," ");	        /* 定数のトークンを読み込む */
      b[i]=atof(input->string);	       	        /* 右辺の定数bに代入する */
      printf("%6.2f\n",b[i]);
   /* if(b[i]==(0.0))		   エラーチェック(0.0)=EROOR  
	{
	  printf("EROOR!!定数データが不正です\n");
	  assert(b[i]!=0.0);	   プログラムを終了する 
         } */
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			read_one_line 
ファイルから１行データを読み込む
エラーがあったらエラーコードを返して終了する
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void read_one_line(FILE *fp,char *line)
{
  /* ファイルから１行読み込む */
  if((fgets(line,MAXWORD,fp))==NULL)
    {	/* １行読み終わったらエラーチェックをする */
      if(feof(fp)!=0)	
	{/* ファイルの中身がないときはエラーメッセージを出力する */
	  printf("EROOR!!ファイルにデータがありません！！\n");
	  assert(feof(fp)==0);	/* プログラムを終了する */
	}
    }
}
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			read_one_token 
	１つのトークンを構文解析する関数
[10.5X1]や[20.0x2]のようなトークンを読み込み
Ｘの係数および（何番目のＸかを表す）Ｘの添字を返す関数
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int read_one_token(char *string,double *num)
{
  int cnt=0;			        /* カウンタ */
  int x_num;			        /* Ｘの添字 */
  char *buf;			        /* 一時記憶用文字配列ポインタ */
  *num=atof(string);			/* 文字型の係数をダブル型に変換する */
  if(*num==(0.0))			/* エラーチェック(0.0)=EROOR */ 
    {
      printf("EROOR!!係数データが不正です\n");
      assert(*num!=0.0);		/* プログラムを終了する */
    }
  while((string[cnt]!='x')&&(string[cnt]!='X'))
    {	                                /* Ｘを発見するまで繰り返す */
      cnt++;
    }
  buf=&string[cnt+1];			/* bufはＸの添字になる */
  x_num = atoi(buf); 			/* 文字型のbufをint型に変換する */
  return x_num;
}
/*-----------------------------------------------------------------------------
　　　関数名：output_title
	　内容：ファイルにタイトルを出力する関数
	　引数：出力ファイル(fout)
-----------------------------------------------------------------------------*/
void output_title(FILE *fout)
{
  int i,j;
  fprintf(fout,"２段階シンプレックス法\n");
  fprintf(fout,"目的関数\n Z =");
  for(i=1;i<=n;i++) fprintf(fout,"%+7.2fX%d",c[i],i);
  if(max_min == max) fprintf(fout," ---> max");
  else  fprintf(fout," ---> min");
  
  for(i=1;i<=m;i++)
    {
      fprintf(fout,"\n条件式%d",i);
      for(j=1;j<=n;j++)
	{
	  if(j==1)
	    fprintf(fout,"%7.2fX%d",a[i][j],j);
	  else
	    fprintf(fout,"%+7.2fX%d",a[i][j],j);
	}
      fprintf(fout," %c= %+7.2f",sign[i], b[i]);
    }
  fprintf(fout,"\n");
}
