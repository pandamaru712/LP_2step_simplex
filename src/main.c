#include "simplex.h"

int main(){
	/* ファイル */
   FILE *fout;     /* FILEポインタ型の*foutを定義する */
	fout = fopen("result.txt","w");

	/* ファイルからデータの読み込み */
	read_data();
	/* タイトルおよび問題の出力 */
	output_title(fout);
	/* シンプレックス法の計算 */
	calculate_simplex(fout);

	fclose(fout);
	printf("\n結果をresult.txtに出力しました\nプログラムを終了します\n");
	return 0;
}
