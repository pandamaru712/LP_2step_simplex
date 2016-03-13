/******************************************************************************
　　　　　　　　LINEAR PROGRAMMING SIMPLEX METHOD改訂版

目的：一般線形計画問題を２段階シンプレックス法を用いて解くプログラム
作成者：h.ago & y.toyota & m.hisai		作成日：2000 12.8
修正者:m.murao                                  修正日: 2001 10.16	       プログラムの使い方：
　あらかじめ入力データを作成しておく
　プログラムを実行するとデータファイルから目的関数および制約条件を読み込む

   Maximize Z=C*X subject to A*X<=b A*X=b A*X>=b

　制約条件の不等号は’＜＝’か’＞＝’か’＝’である
　シンプレックス法を用いて解いた結果をresult.txtに出力する		      								
******************************************************************************/
#include "simplex.h"

/*-----------------------------------------------------------------------------
　　　関数名：calculate_simplex：シンプレックス法の計算
	　内容：シンプレックス演算を行う核となる関数
	　引数：出力ファイル fout
-----------------------------------------------------------------------------*/
void calculate_simplex(FILE *fout)
{
        int rept;                       /* 繰り返し変数 */
	int pivot_retu;	                /*ピボットの列番号 */
	int pivot_gyo;	                /*ピボットの行番号 */
	double pivot;                   /*ピボットの値 */
	double theta[CONSTRAIN_MAX];	/* シータの値 */
	double omega[CONSTRAIN_MAX];	/* オメガ */
	init();		         /* 基底形式の作成 */
	simplex_criterion();     /* シンプレックス基準の初期設定 */
	/* 第１フェーズの計算 */
	rept = 0;     
	reset_omega(omega);     /* ωの初期設定 */
	while(get_first(omega) == yes)	/*ωに負の値がある間*/
	{
		/*ωの絶対値が最大となる列番号を求める*/
		pivot_retu = get_max_omega(omega);
		/*θを求める*/
		calculate_theta(theta,pivot_retu);
		pivot_gyo = get_min_theta(theta);     /* θが最小となる行の番号を求める */
		/* 出力 */
		output_tableau(fout,theta);
		output_omega(fout,omega);
		/* 基底変換 */
		base_change(pivot_gyo,pivot_retu,&pivot);  /* &はアドレスを表す */
		/*１サイクル進める*/
		pivoting(pivot_gyo,pivot_retu,pivot,omega);
		/* whileループの強制終了 */
		if(++rept > nn) break; 
	}
	/* 第２フェーズの計算 */
	rept=0;
	while(get_second() == yes)
	{
	        rept+=1;
		/* c の絶対値が最大となる列番号を求める */
		pivot_retu = get_max_c();
		/* θを求める */
		calculate_theta(theta,pivot_retu);
		pivot_gyo = get_min_theta(theta);     /* θが最小となる行の番号を求める */
		/* 出力 */
		output_tableau(fout,theta);
		if(rept==1) output_omega(fout,omega);
		/* 基底変換 */
		base_change(pivot_gyo,pivot_retu,&pivot);
		/*１サイクル進める */
		pivoting(pivot_gyo,pivot_retu,pivot,omega);
		/* whileループの強制終了 */
		if(rept > nn) break;
	}
	/* 最後θの値は表示されないように負の値にする */
	reset_theta(theta);
	output_tableau(fout,theta);
	if( rept == 0) output_omega(fout,omega);
	if( max_min == min ) z *= -1;
	
	/* 結果の出力 */
	output_result(fout);	
}
/*-----------------------------------------------------------------------------
　　　関数名：init
	　内容：基底形式の作成
	　引数：なし
-----------------------------------------------------------------------------*/
void init()
{
	int i,j;
	nn = n;    
	/* b<0の場合の前処理 */
	for(i=1; i<=m; i++)
	  {
            if(b[i]<0.0)             
	       {
                 for(j=1; j<=n; j++)
		    a[i][j] = -a[i][j];
                 b[i]=-b[i];
		 if(sign[i]=='>') sign[i]='<';
		 else if(sign[i]=='<') sign[i]='>';
	       }
	    for(j=n+1; j<=n+m; j++) a[i][j]=0.0;
	  }	
	/* 余裕変数の追加 */
	for(i=1;i<=m;i++)  
          if(sign[i]=='>')
		{
				nn+=1;
				a[i][nn]=-1.0;
				sign[i]='=';
				c[nn]=0.0;
		}
	/* スラック変数の追加 */
	artivari_num = 0;
	for(i=1; i<=m; i++){
		for(j=nn+1; j<=nn+m; j++)
			a[i][j]=0.0;
		a[i][nn+i]=1.0;
		c[nn+i]=0.0;
		basic_variable[i]= nn+ i;

		if(sign[i]=='=')
		{
			/* 人為変数の数とその番号 */
			artivari_num += 1;
			artivari[artivari_num] = nn+ i;
		}
	}
	nn = nn+m;
}
/*-----------------------------------------------------------------------------
　　　関数名：simplex_criterion
	　内容：シンプレックス基準の初期設定を行う関数
	　引数：なし
-----------------------------------------------------------------------------*/
void simplex_criterion()
{
	int j;
	if(max_min == max){
		for(j=1; j<=nn; j++)
			c[j] *= -1.0;
	}
	z = 0.0;
	/* minであれば「z'」で出力する */
}
/*-----------------------------------------------------------------------------
　　　関数名：reset_omega
	　内容：ωの初期設定を行う関数
	　引数：omega
-----------------------------------------------------------------------------*/
void reset_omega(double omega[])
{
	int i,j;
	for(j=0; j<=nn; j++)
		omega[j] = 0.0;
	for(i=1; i<=m; i++)
	{
		if(sign[i] == '=')
		{
			omega[0] -= b[i];
			for(j=1; j<=nn - m; j++)
				omega[j] -= a[i][j];	
		}
	}
}
/*-----------------------------------------------------------------------------
　　　関数名：get_first
	　内容：第１フェーズの計算を実行するかどうかを判定する関数
	　引数：omega
-----------------------------------------------------------------------------*/
int get_first(double omega[])
{
	int j;
	/*ωにマイナスがあるかチェック*/
	for(j=1; j<=nn; j++)
		if(omega[j] < -ZI) return yes;
	/*ωに値がなければ第２フェーズへ行く*/
	return no;
}
/*-----------------------------------------------------------------------------
　　　関数名：get_max_omega
	　内容：ωの絶対値が最大となる列番号を求めそれを戻す関数
	　引数：omega
-----------------------------------------------------------------------------*/
int get_max_omega(double omega[])
{
	int j,abs_max;
	abs_max = 1;
	for(j=2; j<=nn; j++)
		if( omega[abs_max] > omega[j])
			abs_max = j;
	return abs_max;
}
/*-----------------------------------------------------------------------------
　　　関数名：calculate_theta
	　内容：θを求める関数
	　引数：theta,pivot_retu
-----------------------------------------------------------------------------*/
void calculate_theta(double theta[],int pivot_retu)
{
	int i;
	theta[0] = INF;
	for(i=1; i<=m; i++){
		if(a[i][pivot_retu] > 0.0)
			theta[i] = b[i] / a[i][pivot_retu];
		else theta[i]=-1.0;
	}
}
/*-----------------------------------------------------------------------------
　　　関数名：get_min_theta
	　内容：θが最小となる行の番号を求めそれを戻す関数
	　引数：theta
-----------------------------------------------------------------------------*/
int get_min_theta(double theta[])
{
	int i,minimal;
	minimal = 0;
	for(i=1; i<=m; i++)
	{
		if( theta[i] < 0.0 )	continue;
		if( theta[minimal] > theta[i] )
			minimal = i;
	}
	return minimal;
}
/*-----------------------------------------------------------------------------
　　　関数名：base_change
	　内容：基底変換を行う関数
	　引数：pivot_gyo,pivot_retu,pivot
-----------------------------------------------------------------------------*/
void base_change(int pivot_gyo,int pivot_retu,double *pivot)
{
	*pivot = a[pivot_gyo][pivot_retu];
	basic_variable[pivot_gyo] = pivot_retu;
}	
/*-----------------------------------------------------------------------------
　　　関数名：pivoting
	　内容：１サイクル進める関数
	　引数：pivot_gyo,pivot_retu,pivot,omega
-----------------------------------------------------------------------------*/
void pivoting(int pivot_gyo,int pivot_retu,double pivot,double omega[])
{
	/* ピボットの存在する行の演算 */
	int i,j,gyo,retu;
	double multi;
	for(j=1; j<=nn; j++)
		a[pivot_gyo][j] /= pivot;
	b[pivot_gyo] /= pivot;
	/* ピボットの存在しない行の演算 */
        gyo = pivot_gyo;
	retu = pivot_retu;
	for(i=1; i<=m; i++){
		if(i == gyo ) continue;
		/* 各行への倍率を求める */
		multi = a[i][retu];
		for(j=1; j<=nn; j++)
			a[i][j] -= multi * a[gyo][j];
		b[i] -= multi * b[gyo];
	}
	/* 目的関数部分 */
	multi = c[retu];
	for(j=1; j<=nn; j++)
		c[j] -= multi * a[gyo][j];
	z -= multi * b[gyo];
	multi = omega[retu];
	omega[0] -= multi * b[gyo];
	for(j=1;j<=nn;j++)
		omega[j] -= multi * a[gyo][j];
}
/*----------------------------------------------------------------------------
　　　関数名：get_second
	　内容：第２フェーズの計算を実行するかどうかを判定する関数
	　引数：なし
-----------------------------------------------------------------------------*/
int get_second(){
	int j,k;
	k=1;
	for(j=1;j <= nn;j++)
	{
		if(j == artivari[k])
		{
			k+=1;
			if(k > artivari_num) k-=1;
			continue;
		}
		else if(c[j] < -ZI) return yes;
	}
	return no;
}
/*-----------------------------------------------------------------------------
　　　関数名：get_max_c
	　内容：cの絶対値が最大となる列番号を求めそれを戻す関数
	　引数：なし
-----------------------------------------------------------------------------*/
int get_max_c(){
	int j,k,abs_max;
	abs_max = 1;
	k=1;
	for(j=2; j <= nn ;j++)
	{
		if(j == artivari[k])
		{
			k+=1;
			if(k > artivari_num) k-=1;
			continue;
		}
		else if(c[abs_max] > c[j]) abs_max=j;
	}
	return abs_max;
}
/*-----------------------------------------------------------------------------
　　　関数名：reset_theta
	　内容：θを負の値にする関数
	　引数：theta
-----------------------------------------------------------------------------*/
void reset_theta(double theta[]){
	int j;
	for(j=1; j<=m; j++)
		theta[j] = -1.0;
}
/******************************************************************************
----------------- 出力 --------------------------------------------------------
*******************************************************************************
　　　関数名：output_tableau
	　内容：１サイクル分のタブローを出力する関数
	　引数：出力ファイル fout,theta
-----------------------------------------------------------------------------*/
void output_tableau(FILE *fout,double theta[])
{
  static int count = 0;
  int i,j,k;
  fprintf(fout,"\n///// サイクル = %d /////\n",count);
  
  /* タブローの見出しの出力 */
  fprintf(fout,"基底変数とその値  ");
  for(j=1; j<=n; j++)
    fprintf(fout,"  Ｘ%d    ",j);
  for(j=1; j<=nn - n; j++)
    fprintf(fout,"  λ%d    ",j);
  fprintf(fout,"   θ\n");
  fprintf(fout,"  -----------------------");
  for(j=1; j<=nn; j++)
    {
      fprintf(fout,"---------");
    }
  fprintf(fout,"\n");
  
  /* 制約条件部の出力 */
  for(i=1; i<=m; i++)
    {
      if(basic_variable[i] <= n) fprintf(fout,"  Ｘ%d",basic_variable[i]);
      else fprintf(fout,"  λ%d",basic_variable[i] -n);
      fprintf(fout," %9.2f",b[i]);
      for(j=1; j<=nn; j++)
	fprintf(fout,"%9.2f",a[i][j]);
      if(theta[i] >= 0.0)
	fprintf(fout,"%10.2f\n",theta[i]);
      else	fprintf(fout,"\n");
    }
  
  /* 目的関数部の出力 */
  fprintf(fout,"  -----------------------");
  for(j=1; j<=nn; j++)
    fprintf(fout,"---------");
  fprintf(fout,"\n");
  if(max_min == max) fprintf(fout,"  Ｚ ");
  else fprintf(fout,"  Z' ");
  fprintf(fout," %9.2f",z);
  k=1;
  for(j=1; j<=nn; j++)
    {
      if(j != artivari[k])
	fprintf(fout,"%9.2f",c[j]);
      else
	{
	  fprintf(fout,"         ");
	  k+=1;
	  if(k > artivari_num) k-=1;
	}
    }
  fprintf(fout,"\n");
  
  ++count;
}
/*-----------------------------------------------------------------------------
　　　関数名：output_omega
	　内容：ωを出力する関数
	　引数：出力ファイル fout,omega
-----------------------------------------------------------------------------*/
void output_omega(FILE *fout,double omega[])
{
  int j;
  fprintf(fout,"  ω  ");
  for(j=0; j<=nn; j++)
    fprintf(fout,"%9.2f",omega[j]);
  fprintf(fout,"\n");
}

/********** 結果(result) ******************************************************
　　　関数名：output_result
	　内容：解を出力する関数
	　引数：出力ファイル(fout)
-----------------------------------------------------------------------------*/
void output_result(FILE *fout)
{
  int i;
  fprintf(fout,"\n/////////// 解 //////////////\n");
  fprintf(fout,"Z   = %f\n",z);
  for(i=1; i<=n; i++)
    fprintf(fout,"X%d  = %f\n",i ,result_read(i));
  for(i=n+1; i<=nn; i++)
    fprintf(fout,"λ%d = %f\n",i-n ,result_read(i));
}
/*-----------------------------------------------------------------------------
　　　関数名：resulut_read
	　内容：出力したい変数が基底変数に含まれているかチェックし，
			含まれていればその値を，そうでなければゼロを返す関数
	　引数：出力する変数の番号
-----------------------------------------------------------------------------*/
double result_read(int i)
{
  int j;
  for(j=1; j<=m; j++)
    if(i==basic_variable[j] )
      return b[j];
  return 0.0;
}
