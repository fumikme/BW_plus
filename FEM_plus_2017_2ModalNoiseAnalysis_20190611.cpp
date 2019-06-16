/***  更新履歴  2019/05/04  *****///							ここから印刷切れます↓
/************************************************************************************/
/*-*		     【  -*-*-*-*-*-    関心毎の分離    -*-*-*-*-*-  】				  *-*/
/*-*																			  *-*/
/*-*			モジュール化してそれぞれ10行くらいにしておくと，			      *-*/
/*-*	  実装する時にも分離した小さい単位のロジックに思考を集中できるので、	  *-*/
/*-*	                  コードが書きやすくなります。						      *-*/
/*-*																			  *-*/
/************************************************************************************/
/* iostreamが読み込めない際は，Windows SDKバージョンを更新してください */
#define _USE_MATH_DEFINES
#define NR_END 1
#define FREE_ARG char*
#pragma warning(disable:4996)

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <direct.h>
#include <sys\stat.h>
using namespace std;

/* 諸定数の定義 */
#define PI 3.141592653	
#define C 2.99792458e8	
#define Epsilon0 8.8541878e-12 
#define Mu0 PI*4.0e-7	

/* 出力副プログラムの読み込み */
void sub_routine_4_2_1_mkdir(char dir[], int mrow);
void sub_routine_4_2_2(char dir[], int mrow, int noffV2F, int noffF2F, double ****TotalCEF,double ***MPDoutput1st, double ***MPDoutput2nd);
void sub_routine_4_2_3(char dir[], int mrow, int noffV2F, int noffF2F, double ****TotalCEF);
void sub_routine_4_2_4(char dir[], int mrow, int noffV2F, int noffF2F, double ****TotalCEF);

/* 関数副プログラムの読み込み */
// 屈折率波長微分
double dndl(double lamda, double n_lamda, int mater);
// 係数行列要素計算関数
void S_matrix(double *a, double *b, double *q, int m, int n, double v, double w, double D);
// 改訂コレスキー分解および改訂コレスキー分解による求解
void mcholesky(double *a, double *b, double *ML, double *MD, int m, int n);
void mcholesky_sol(double *ML, double *MD, double *R, int m, int n);
// 逆べき乗法の初期ベクトル計算・解ベクトル規格化
void R0(double *MD, double *R, int m, int n);
void R_norm(double *R, int n);
// 固有値計算
double Eigen(double *R, double *a, double *b, int n, int m);
//  群遅延計算用関数
double dbdk_bunbo(double *R, double D, double w, int m, int n);
double dbdk_bunshi(double *R, double *qg, double D, double w, int m, int n);
// 畳み込み積分
double *convolution(int xmax, double *P, int LT, double *pulse);
//																 ここから印刷切れます↓
/* 数学関数・多次元配列の読み込み */
// 第1種変形Bessel関数 In(x)・第2種変形Bessel関数 Kn(x)
double bessi0(double x), bessi1(double x);
double bessk0(double x), bessk1(double x), bessk(int n, double x);
// 1次元配列（整数型）
int  *dintvector(int i, int j);
void free_dintvector(int *a, int i);
void init_intvector(int *a, int nr1, int nr2);
// 1次元配列（実数型）
double *drealvector(int i, int j);
void free_drealvector(double *a, int i);
void init_realvector(double *a, int nr1, int nr2);
// 2次元配列（float型）
float **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(float **a, long nrl, long nrh, long ncl, long nch);
void init_matrix(float **a, int nr1, int nr2, int nl1, int nl2);
// 2次元配列（double型）
double  **dmatrix(int nr1, int nr2, int nl1, int nl2);
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
void init_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
// 3次元配列
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***a, long nrl, long nrh, long ncl, long nch,	long ndl, long ndh);
void init_d3tensor(double ***a, int nr1, int nr2, int nl1, int nl2, int np1, int np2);
// 4次元配列
double ****d4tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long ntl, long nth);
void free_d4tensor(double ****a, long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long ntl, long nth);
void init_d4tensor(double ****a, long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long ntl, long nth);

/* d3tensor の第2、第3成分の初期化*/
void init_d3tensor_nfnz(double ***a, int nr, int nl1, int nl2, int np1, int np2);

/* エラー出力関数の読み込み */ 
void nrerror(char error_text[]);
/* エラー出力関数の読み込み */
void mkdir(char dirname[]);

/* パラメータの定義 
*  m: モード次数( TE&TM ~ 1, EH ~ n+1, HE ~ n-1), l: 動径方向モード次数, n: 方位角モード次数
*  mater: 材料ID, profile: 屈折率分布入力方法ID (0: べき乗分布, 1: 測定分布入力)
*  NLP: LPモード数, NLP0: 0次方位角モードLPモード数, Ntotal: 全モード数
*  N: 動径座標rのコア内空間ステップ数, Nclad: 解析領域空間ステップ数, Nbeta: 伝搬定数の分割数
*  Nwkb: WKB近似により求めたLPモード数, Ntmin: 最小群遅延時間ステップ数, Ntmax: 最大群遅延時間ステップ数
*  lamda: 波長, k: 波数, omega: 角周波数, delta: 比屈折率, NA: 開口数, aa: 動径座標規格化サイズ
*  A: コア半径, AA: 解析領域半径, g: べき乗分布指数, n1: コア中心の屈折率, n0: クラッドの屈折率
*  L: ファイバ長, Tv: インパルス応答時間軸刻み幅, dr: 動径座標刻み幅,
*  v: 規格化周波数, w: 規格化伝搬定数, D: 規格化コア径
*  tau: 群速度, beta: 伝搬定数β, dbeta: 伝搬定数刻み幅
*  a: 係数行列Sの副対角要素格納配列, b: 係数行列Sの対角要素格納配列
*  GI: コア内屈折率, q: 規格化コア内屈折率, qg: コア内屈折率分散パラメータ
*  R: 規格化横方向電場成分, Rb: 逆べき乗法用入れ子配列, Rinf: 解析領域
*  ML: LDU分解係数行列のL行列副対角要素, MD: LDU分解係数行列のD行列対角要素
*  de, df: 逆べき乗法における収束判別パラメータ, eig: 固有値
*  taumin: 最小単位群遅延, taumax: 最大単位群遅延
*  fmin: 評価最小周波数, fmax: 評価最大周波数, dfrq: 評価周波数分解能
*  Hw: 伝達関数要素入れ子変数 (ReHw: 実部, ImHw: 虚部), bw: -3dB帯域幅
*  GI: 屈折率分布, q: 屈折率波長微分, qg: n*d(kn)/dk
*  Rlp: 全LPモードのモードフィールド, Mbeta: 各モードの伝搬定数
*  r0: オフセット（励振位置），w0: ガウシアンビーム半径，dx: 重なり積分の空間分解能, Nxy: 解析領域の分割数 */

//																ここから印刷切れます↓
int main(void)			 // tensor=(page, column, row)
{	
// 0.1  パラメータの設定
#define Nvtotal 6
#define Nftotal 80		 // #define Nftotal 79  #define Nftotal 90

#define OFFres_V2F	     1
#define OFFrange_V2F     20
#define OFFres_F2F	     1
#define OFFrange_F2F     20

#define fiberlength      30  	// unit: [m]
#define outputlength     5		//  〃 : [m]
#define	calcultlength    1 		//  〃 : [m]    …【 1 に固定】
#define z_res            0.001	//  〃 : [m]
#define fiber		    "XX3"	//【Silica2】【GigP】【XX3】
#define dirnamecoeff    "Coupling_coefficient"
#define dirMPD_V2F		"MPD_V2F"
#define dirMPD_F2F		"MPD_F2F"
#define filenameV2F		"MSB_V2F_Offset_"
#define filenameF2F		"TotalCEF_F2F_Offset_"
#define unit            "um"
#define directoryMPDout "CEFandMPDperVCSELMS(zdep)"
#define dirTotalCEF     "AllofTotalCEF"

#define dirdebug        "dirdebug"

	// 0.1  変数の宣言と初期化
	FILE *fhuv, *fMPDinput, *fMPD_F2Fcoupling, *fsumhuv, *fallCEF;
	FILE *fdebug, *fdebug1, *fdebug2;

	char frpath[256], fibertypecsv[256], trash[65536];
		 frpath[0] =  fibertypecsv[0] =  trash[0] =  '\0';
	char MPDinfilename[256], MPDoutdirectory[256], MPD_F2Fcoupling[256];
		 MPDinfilename[0] =  MPDoutdirectory[0] =  MPD_F2Fcoupling[0] =  '\0';
	char AllofCEFfilename[256], sumhuvcsv[256], fnamedebug[256], filename[256];
		 AllofCEFfilename[0] =  sumhuvcsv[0] =  fnamedebug[0] =  filename[0] =  '\0';
	int  i,j;	int nv,nf,nz;	int noffV2F,noffF2F;	int m1,m2;		int nf_out,nf_in;
		 i=j=0;		nv=nf=nz=0;		noffV2F=noffF2F=0;		m1=m2=0;		nf_out=nf_in=0;
	int  cnt=1;	int nz_cal=0;	int Discriminant=0;
	double *sum_huv, **MPDzcal_bound;
	double **huv;				
	double ***MPDf2f_tensor, ***MPDinput, ***MPDzdep, 
		   ***MPDoutput1st, ***MPDoutput2nd;
	double ****TotalCEF;    //double ***TotalCEF;	 
	
	// 0.2  配列の宣言と初期化
	huv     = dmatrix(0, Nftotal, 0, Nftotal);	
	sum_huv = drealvector(0, Nftotal);		
	MPDf2f_tensor = d3tensor(0, OFFrange_F2F / OFFres_F2F + 1, 0, Nftotal, 0, Nftotal);
	MPDinput = d3tensor(0, Nvtotal, 0, Nftotal, 0, 1);			
	MPDzdep = d3tensor(0, Nvtotal, 0, Nftotal, 0, (int)(calcultlength / z_res));
	MPDzcal_bound = dmatrix(0, Nftotal, 0, (int)(calcultlength / z_res));
	MPDoutput1st = d3tensor(0, Nvtotal, 0, Nftotal, 0, fiberlength / calcultlength);
	MPDoutput2nd = d3tensor(0, Nvtotal, 0, Nftotal, 0, fiberlength / calcultlength);
	TotalCEF = d4tensor(0, Nvtotal, 0, OFFrange_V2F / OFFres_V2F + 1, 0, OFFrange_F2F / OFFres_F2F + 1, 0, fiberlength/calcultlength);
	   init_dmatrix(huv, 0, Nftotal, 0, Nftotal);
	   init_realvector(sum_huv, 0, Nftotal);
	   init_d3tensor(MPDf2f_tensor, 0, OFFrange_F2F / OFFres_F2F + 1, 0, Nftotal, 0, Nftotal);
	   init_d3tensor(MPDinput, 0, Nvtotal, 0, Nftotal, 0, 1);
	   init_d3tensor(MPDzdep, 0, Nvtotal, 0, Nftotal, 0, (int)(calcultlength / z_res));
	   init_dmatrix(MPDzcal_bound, 0, Nftotal, 0, (int)(calcultlength / z_res));
	   init_d3tensor(MPDoutput1st, 0, Nvtotal, 0, Nftotal, 0, fiberlength / calcultlength);
	   init_d3tensor(MPDoutput2nd, 0, Nvtotal, 0, Nftotal, 0, fiberlength / calcultlength);
	   init_d4tensor(TotalCEF, 0, Nvtotal, 0, OFFrange_V2F / OFFres_V2F + 1, 0, OFFrange_F2F / OFFres_F2F + 1, 0, fiberlength / calcultlength);
	
	   // 4.0.2  デバッグ用ディレクトリの作成
	//mkdir(dirdebug);
	sprintf(fnamedebug, "%s\\%s", dirdebug, fiber);				mkdir(fnamedebug);
	sprintf(fnamedebug,"%s\\%s\\MPD2d_F2F", dirdebug, fiber);	mkdir(fnamedebug);
	sprintf(fnamedebug,"%s\\%s\\MPDinput", dirdebug, fiber);	mkdir(fnamedebug);
	sprintf(fnamedebug,"%s\\%s\\MPDzdep", dirdebug, fiber);		mkdir(fnamedebug);
	sprintf(fnamedebug,"%s\\%s\\MPDoutput1st", dirdebug, fiber);mkdir(fnamedebug);
	sprintf(fnamedebug,"%s\\%s\\MPDoutput2nd", dirdebug, fiber);mkdir(fnamedebug);
	sprintf(fnamedebug,"%s\\%s\\TotalCEF", dirdebug, fiber);	mkdir(fnamedebug);

//																ここから印刷切れます↓
	// 1.  結合係数huvを2D配列に格納
	sprintf(fibertypecsv, "%s\\%s.csv", dirnamecoeff, fiber);
	int    tr1 = 0;	 //int ti[20];
	double tr2 = 0.0;	//double td[20];
	if (fhuv = fopen(fibertypecsv, "r")) {
		fscanf(fhuv, "%[^\n]\n", trash);
		fscanf(fhuv, "%[^\n]", trash);
		for (i = 0; i <= Nftotal; i++) {
			for (j = 0; j < i; j++) {
				fscanf(fhuv, "%d,%d,%d,%d,%lf,%d,%d,%d,%d,%lf,%d,%lf,%lf,%lf\n",
					&tr1, &tr1, &tr1, &tr1, &tr2, &tr1, &tr1, &tr1, &tr1, &tr2, &tr1, &tr2, &tr2, &huv[i][j]);
				/*fscanf(fhuv, "%d,%d,%d,%d,%lf,%d,%d,%d,%d,%lf,%d,%lf,%lf,%lf\n",
					&ti[0], &ti[1], &ti[2], &ti[3], &td[4], &ti[5], &ti[6], &ti[7], &ti[8], &td[9], &ti[9], &td[10], &td[11], &huv[i][j]);
				printf("%d, %d, %d, %d, %lf, %d, %d, %d, %d, %lf, %d, %lf, %lf, %lf\t",
					ti[0], ti[1], ti[2], ti[3], td[4], ti[5], ti[6], ti[7], ti[8], td[9], ti[9], td[10], td[11], huv[i][j]);
				printf("%lf\n", huv[i][j]);	*/			}}
		fclose(fhuv); }
	if(fiber=="GigP" || fiber=="XX3") { goto next; }
	if (fhuv = fopen(fibertypecsv, "r")) {
		fscanf(fhuv, "%[^\n]\n", trash);
		fscanf(fhuv, "%[^\n]", trash);
		for (i = 0; i <= Nftotal; i++) {
			for (j = 0; j < i; j++) {
				fscanf(fhuv, "%d,%d,%d,%d,%lf,%d,%d,%d,%d,%lf,%d,%lf,%lf\n",
					&tr1, &tr1, &tr1, &tr1, &tr2, &tr1, &tr1, &tr1, &tr1, &tr2, &tr1, &tr2, &huv[i][j]);	
				//printf("%lf\n", huv[i][j]);		
			}}
		fclose(fhuv); }
	next:
	for (i = 0; i <= Nftotal; i++) {
		for (j = 0; j < i; j++) { huv[j][i] = huv[i][j];				}}
	sprintf(sumhuvcsv, "%s\\sumhuv_%s.csv", dirnamecoeff, fiber);
	if (fsumhuv = fopen(sumhuvcsv, "w")){
		fprintf(fsumhuv, "mode number, sum_huv[number]\n");
		for (i = 0; i <= Nftotal; i++) {
			for (j = 0; j <= Nftotal; j++) { 
				sum_huv[i] += huv[i][j];					}
			fprintf(fsumhuv, "%d, %e\n", i, sum_huv[i]); }
		fclose(fsumhuv);  }
	sprintf(sumhuvcsv, "%s\\huvmatrix_%s.csv", dirnamecoeff, fiber);
	if (fsumhuv = fopen(sumhuvcsv, "w")){
		fprintf(fsumhuv, "mode number, ");
		for (i = 0; i <= Nftotal; i++) { fprintf(fsumhuv, "%d,", i); } fprintf(fsumhuv, "\n");
		for (i = 0; i <= Nftotal; i++) {
			fprintf(fsumhuv, "%d,", i);
			for (j = 0; j <= Nftotal; j++) { 
				fprintf(fsumhuv, "%e,", huv[i][j]);			}
			fprintf(fsumhuv, "\n");}
		fclose(fsumhuv);  }
	printf("1.  結合係数huvを2D配列に格納しました。\n");

	// 2.  F2FのMSBを3D配列に格納
	for (noffF2F = 0; noffF2F < OFFrange_F2F / OFFres_F2F + 1;noffF2F++) {
		sprintf(MPD_F2Fcoupling, "%s\\MPD2d_F2F_%dum_Offset.csv", dirMPD_F2F, noffF2F*OFFres_F2F);	
		if (fMPD_F2Fcoupling = fopen(MPD_F2Fcoupling, "r")) {
			sprintf(fnamedebug, "%s\\%s\\MPD2d_F2F\\MPDf2f_tensor_x2_%dum_Offset.csv", dirdebug, fiber, noffF2F*OFFres_F2F);
			fdebug = fopen(fnamedebug, "w");
			fscanf(fMPD_F2Fcoupling, "%[^\n]\n", trash);
			fprintf(fdebug, "-,");
			for (nf_in = 0; nf_in < Nftotal; nf_in++) { fprintf(fdebug, "%d,", nf_in); }fprintf(fdebug, "\n");
			for (nf_out = 0; nf_out < Nftotal; nf_out++) {
				fscanf(fMPD_F2Fcoupling, "%lf, ", &tr2);
				fprintf(fdebug, "%d,", nf_out);
				for (nf_in = 0; nf_in < Nftotal; nf_in++) {
					fscanf(fMPD_F2Fcoupling, "%lf, ", &MPDf2f_tensor[noffF2F][nf_out][nf_in]);	
					MPDf2f_tensor[noffF2F][nf_out][nf_in] = MPDf2f_tensor[noffF2F][nf_out][nf_in]/100.0;
					fprintf(fdebug, "%lf,", MPDf2f_tensor[noffF2F][nf_out][nf_in]);	}
				fprintf(fdebug, "\n");		}fclose(fdebug);
			}}
	printf("2.  F2FのMSBを3D配列に格納しました。\n");

	// 3.  MPD出力用ディレクトリを作成
	struct stat statBuf;
	char dir[256];
	char xdir[256];
	sprintf(dir, "%s_%s", fiber, directoryMPDout);
	if (stat(dir, &statBuf) != 0){ 
		if (_mkdir(dir) != 0) {
			  printf("ディレクトリ %s の作成に失敗しました。\n", dir);
			  system("pause"); }
		else{ printf("3.  MPD出力用ディレクトリ %s の作成に成功しました。\n", dir);	}}
	else{ printf("3.  MPD出力用として、既に作成されたディレクトリ %s を使用します。\n", dir); }
	sprintf(xdir, "%s\\%s_x", dir, dirTotalCEF);	
	mkdir(xdir);

//																ここから印刷切れます↓
	// 4.  MSBとモード結合を考慮して，3D配列MPDの変化を追跡
	// 4.0.1  出力ファイルのヘッダ作成（outputlength毎に加え、1mに関しても作成）
	char fname[256];
	for (int m1 = 0; m1 <= fiberlength / outputlength; m1++) {
		sprintf(fname, "%s\\%s_z_%dm.csv", dir, dirTotalCEF, m1*outputlength);
		fallCEF = fopen(fname, "w");
		fprintf(fallCEF, "x1, x2, ");
		for (nv = 0; nv < Nvtotal; nv++) { 
			fprintf(fallCEF, "MS%d, ", Nvtotal - nv); }
		fprintf(fallCEF, "\n");
		fclose(fallCEF);		}
	if(outputlength>1){
		sprintf(fname, "%s\\%s_z_%dm.csv", dir, dirTotalCEF, 1);
		fallCEF = fopen(fname, "w");
		fprintf(fallCEF, "x1, x2, ");
		for (nv = 0; nv < Nvtotal; nv++) { 
			fprintf(fallCEF, "MS%d, ", Nvtotal - nv); }
		fprintf(fallCEF, "\n");
		fclose(fallCEF);		
	}
	for (noffV2F = 0; noffV2F < OFFrange_V2F / OFFres_V2F + 1; noffV2F++) {
		for (noffF2F = 0; noffF2F< OFFrange_F2F / OFFres_F2F + 1; noffF2F++) {
			sprintf(fname, "%s\\%s_x\\%s_x1=%dum_x2=%dum.csv", dir, dirTotalCEF, dirTotalCEF, noffV2F*OFFres_V2F, noffF2F*OFFres_F2F);
			fallCEF = fopen(fname, "w");
			fprintf(fallCEF, "x1, x2, z, ");
			for (nv = 0; nv < Nvtotal; nv++) {
				fprintf(fallCEF, "MS%d, ", Nvtotal - nv);
			}
			fprintf(fallCEF, "\n");
			fclose(fallCEF);		}}

	// 4.0.2  レーザと光ファイバ間結合におけるMSBの読み込み
	for (noffV2F = 0; noffV2F < OFFrange_V2F / OFFres_V2F + 1; noffV2F++){
		sprintf(MPDinfilename, "%s\\MPD2d_V2F_%dum_Offset.csv", dirMPD_V2F, noffV2F*OFFres_V2F);
		sprintf(fnamedebug, "%s\\%s\\MPDinput\\MPDinput_x1_%dum_Offset.csv", dirdebug, fiber, noffV2F*OFFres_V2F);
		fdebug = fopen(fnamedebug, "w");
		if (fMPDinput = fopen(MPDinfilename, "r")) {
			fscanf(fMPDinput, "%[^\n]\n", trash);
			fprintf(fdebug, "-,");
			for (nv = 0; nv < Nvtotal; nv++) { fprintf(fdebug, "%d,", nv); }fprintf(fdebug, "\n");
			for (nf = 0; nf <= Nftotal; nf++) {
				fscanf(fMPDinput, "%d, %d, %d, ", &tr1, &tr1, &tr1);
				fprintf(fdebug, "%d,", nf);
				for (nv = 0; nv < Nvtotal; nv++) {
					fscanf(fMPDinput, "%lf,", &MPDinput[nv][nf][0]);	
					MPDinput[nv][nf][0] = MPDinput[nv][nf][0] /100.0;
					fprintf(fdebug, "%lf,", MPDinput[nv][nf][0]);
				}fprintf(fdebug, "\n");}fclose(fdebug);}
		
	// 4.1  スプリットステップ法により，モード結合によるMPD変化を計算
		for(nv=0;nv<Nvtotal;nv++){
			for (nz_cal = 0; nz_cal < fiberlength / calcultlength; nz_cal++) {
				init_d3tensor_nfnz(MPDzdep, nv, 0, Nftotal, 0, (int)(calcultlength / z_res));
				for (nz = 0;nz <= calcultlength / z_res;nz++) {
					for (nf = 0;nf < Nftotal;nf++) {
						if (nz == 0 && nz_cal == 0) { MPDzdep[nv][nf][nz] = MPDinput[nv][nf][0]; 	
													  MPDoutput1st[nv][nf][0] = MPDinput[nv][nf][0];	}
						else if (nz == 0) { MPDzdep[nv][nf][nz] = MPDoutput1st[nv][nf][nz_cal]; }
						else { 
							for (int nf_zminus = 0; nf_zminus <= Nftotal; nf_zminus++) {
								if (nf != nf_zminus) {
									MPDzdep[nv][nf][nz] += MPDzdep[nv][nf_zminus][nz - 1] * huv[nf][nf_zminus] * z_res;	
								}
								else if (nf == nf_zminus) { 
									MPDzdep[nv][nf][nz] += MPDzdep[nv][nf_zminus][nz - 1] * (1.0 - sum_huv[nf_zminus] * z_res);	
								}}}
						//if (nf == 0) { printf("%lf\t", MPDzdep[nv][nf][nz]);}	//printf
				}}									
				//printf("\n z=%f\n", outputlength*m1);
				for (nf = 0;nf <= Nftotal;nf++) {			
					MPDoutput1st[nv][nf][nz_cal+1] = MPDzdep[nv][nf][(int)(calcultlength / z_res)];
					//if (m2 == 0) { printf("%d\t%lf\t", m1*outputlength, MPDoutput1st[nv][m2][m1]); }
				}
				if (nz == (int)(1 / z_res) && nz % (int)(outputlength / z_res) != 0) { Discriminant = 1; }}

			sprintf(fnamedebug, "%s\\%s\\MPDzdep\\MPDzdep_mode%d(MS%d)_x1_%dum.csv", dirdebug, fiber, nv+1, Nvtotal-nv, noffV2F*OFFres_V2F);
			fdebug1 = fopen(fnamedebug, "w");
			fprintf(fdebug1, "z_res=%dm,", (int)z_res);
			for (nz_cal = 0; nz_cal <= fiberlength / calcultlength; nz_cal++) { 
				fprintf(fdebug1, "%lf,", (double)nz_cal*calcultlength); }
			fprintf(fdebug1, "\n");
			for (nf = 0; nf < Nftotal; nf++){
				fprintf(fdebug1, "%d,", nf);
				for (nz_cal = 0; nz_cal <= fiberlength / calcultlength; nz_cal++) { 
					fprintf(fdebug1, "%lf,", MPDoutput1st[nv][nf][nz_cal]); }
				fprintf(fdebug1, "\n"); 			}	
			sprintf(fnamedebug, "%s\\%s\\MPDoutput1st\\MPDoutput1st_mode%d(MS%d)_x1_%dum_Offset.csv", dirdebug, fiber, nv+1, Nvtotal-nv, noffV2F*OFFres_V2F);
			fdebug2 = fopen(fnamedebug, "w");
			fprintf(fdebug2, "-,");
			for (nz = 0; nz <= fiberlength / calcultlength; nz++) { fprintf(fdebug2, "%lf,", (double)nz*calcultlength); }
			fprintf(fdebug2, "\n");
			for (nf = 0; nf < Nftotal; nf++){
				fprintf(fdebug2, "%d,", nf);
				for (nz_cal = 0; nz_cal <= fiberlength / calcultlength; nz_cal++) {
					fprintf(fdebug2, "%lf,", MPDoutput1st[nv][nf][nz_cal]); }
				fprintf(fdebug2, "\n"); }	
			for (nf = 0; nf < Nftotal; nf++){
				fprintf(fdebug1, "%d,", nf);
				for (nz_cal = 0; nz_cal <= fiberlength / calcultlength; nz_cal++) { 
					fprintf(fdebug1, "%lf,", MPDoutput1st[nv][nf][nz_cal]); }
				fprintf(fdebug1, "\n"); 			}	
			fclose(fdebug1);
			fclose(fdebug2);
		}
		printf("    4.%d  x1=「%dum」の際の、MPDのz依存性を計算しました。\n", cnt, noffV2F*OFFres_V2F); cnt += 1;
		
	// 4.2  F2FのMSBを用いて，F2Fに軸ズレが生じる際のモード雑音の指標を計算
		for(noffF2F=0; noffF2F< OFFrange_F2F / OFFres_F2F + 1; noffF2F++){
			init_d3tensor(MPDoutput2nd, 0, Nvtotal, 0, Nftotal, 0, fiberlength / calcultlength);
			for (m1 = 0; m1 <= fiberlength / calcultlength; m1++) {
				for(nv=0;nv< Nvtotal;nv++){
					for (nf_out = 0; nf_out < Nftotal; nf_out++) {
						for (nf_in = 0; nf_in < Nftotal; nf_in++) {
							MPDoutput2nd[nv][nf_out][m1] += MPDoutput1st[nv][nf_in][m1] * MPDf2f_tensor[noffF2F][nf_out][nf_in]; }
						//printf("%f\n", MPDoutput2nd[nv][nf_out][m1]);
					}
					for (nf = 0; nf < Nftotal; nf++) { 
						TotalCEF[nv][noffV2F][noffF2F][m1] += MPDoutput2nd[nv][nf][m1]; }
				}}	
			for (nv = 0; nv < Nvtotal; nv++) {
				sprintf(fnamedebug, "%s\\%s\\MPDoutput2nd\\MPDoutput2nd_mode%d(MS%d)_x1_%dum_x2_%dum_Offset.csv", dirdebug, fiber, nv+1, Nvtotal-nv, noffV2F*OFFres_V2F, noffF2F*OFFres_F2F);
				fdebug1 = fopen(fnamedebug, "w");
				fprintf(fdebug1, "-,");
				for (nz = 0; nz <= fiberlength / calcultlength; nz++) { fprintf(fdebug1, "%lf,", (double)nz*calcultlength); }
				fprintf(fdebug1, "\n");
				for (nf = 0; nf < Nftotal; nf++){
					fprintf(fdebug1, "%d,", nf);
					for (nz_cal = 0; nz_cal <= fiberlength / calcultlength; nz_cal++) {
						fprintf(fdebug1, "%lf,", MPDoutput2nd[nv][nf][nz_cal]); }
					fprintf(fdebug1, "\n"); 
				}
				sprintf(fnamedebug, "%s\\%s\\TotalCEF\\TotalCEF_mode%d(MS%d)_x1_%dum_x2_%dum_Offset.csv", dirdebug, fiber, nv+1, Nvtotal-nv, noffV2F*OFFres_V2F, noffF2F*OFFres_F2F);
				fdebug2 = fopen(fnamedebug, "w");
				fprintf(fdebug2, "-,");
				for (nz = 0; nz <= fiberlength / calcultlength; nz++) { fprintf(fdebug2, "%lf,", (double)nz*calcultlength); }
				fprintf(fdebug2, "\n");
				for (nf = 0; nf < Nftotal; nf++){
					fprintf(fdebug2, "%d,", nf);
					for (nz_cal = 0; nz_cal <= fiberlength / calcultlength; nz_cal++) {
						fprintf(fdebug2, "%lf,", MPDoutput2nd[nv][nf][nz_cal]); }
					fprintf(fdebug2, "\n"); }
				fclose(fdebug1);
				fclose(fdebug2);			}
			for (m1 = 0; m1 <= fiberlength / outputlength; m1++) {
				sub_routine_4_2_1_mkdir(dir, m1*outputlength);
				sub_routine_4_2_2(dir, m1*outputlength, noffV2F, noffF2F, TotalCEF, MPDoutput1st, MPDoutput2nd);
				sub_routine_4_2_3(dir, m1*outputlength, noffV2F, noffF2F, TotalCEF);	
				sub_routine_4_2_4(xdir, m1*outputlength, noffV2F, noffF2F, TotalCEF);	}
			if(outputlength >1 || calcultlength<=1) {	
				sub_routine_4_2_1_mkdir(dir, 1);
				sub_routine_4_2_2(dir, (int)(1 / calcultlength), noffV2F, noffF2F, TotalCEF, MPDoutput1st, MPDoutput2nd);
				sub_routine_4_2_3(dir, (int)(1 / calcultlength), noffV2F, noffF2F, TotalCEF);
				sub_routine_4_2_4(xdir, (int)(1 / calcultlength), noffV2F, noffF2F, TotalCEF);		}}
		printf("    4.%d  x1=「%dum」の際の、モード雑音の指標「CEFms」のx2依存性とそのz依存性を計算しました。\n", cnt, noffV2F*OFFres_V2F); cnt += 1;}

	printf("4.  モード雑音の指標「CEFms」の、オフセットx1、x2およびz依存性を計算しました。\n");
	printf("------------------------- Modal_Noise_Analysis_w_wo_MCの終了-------------------------\n\n");
	system("pause");
	return 0;
}

//																ここから印刷切れます↓
// 4.2.1  CEFとMPDを、オフセット量x1，x2の組み合わせ毎に出力
void sub_routine_4_2_1_mkdir(char dir[], int mrow) {
	char zdir[256];	zdir[0] = '\0';
	sprintf(zdir, "%s\\z_%dm", dir, mrow);
	mkdir(zdir);
}
// 4.2.2  CEFとMPDを、オフセット量x1，x2の組み合わせ毎に出力
void sub_routine_4_2_2(char dir[], int mrow, int noffV2F, int noffF2F, 
	double ****TotalCEF, double ***MPDoutput1st, double ***MPDoutput2nd) {
	char str_sub[256];	str_sub[0] = '\0';
	FILE *fMPDoutput;	int nv, nf;	
	sprintf(str_sub, "%s\\z_%dm\\MPDoutput_x1_%dum_x2_%dum.csv", 
		dir, mrow, noffV2F*OFFres_V2F, noffF2F*OFFres_F2F);
	if (fMPDoutput = fopen(str_sub, "w")) {
		fprintf(fMPDoutput, "Modal noise index, ");
		for (nv = 0; nv < Nvtotal; nv++) { fprintf(fMPDoutput, "MS%d, ", Nvtotal - nv); }
		fprintf(fMPDoutput, "\n");	fprintf(fMPDoutput, "TotalCEF, ");
		for (nv = 0; nv < Nvtotal; nv++) { 
			fprintf(fMPDoutput, "%lf, ", TotalCEF[nv][noffV2F][noffF2F][mrow]); }
		fprintf(fMPDoutput, "\n");	fprintf(fMPDoutput, "MPDoutput1st, ");
		for (nv = 0; nv < Nvtotal; nv++) { 
			fprintf(fMPDoutput, "MS%d, ", Nvtotal - nv); }
		fprintf(fMPDoutput, "\n");
		for (nf = 0; nf < Nftotal; nf++) {
			fprintf(fMPDoutput, "%d, ", nf);
			for (nv = 0; nv < Nvtotal; nv++) { 
				fprintf(fMPDoutput, "%lf, ", MPDoutput1st[nv][nf][mrow]); }
			fprintf(fMPDoutput, "\n");		}
		fprintf(fMPDoutput, "\n");	fprintf(fMPDoutput, "MPDoutput2nd, ");
		for (nv = 0; nv < Nvtotal; nv++) { fprintf(fMPDoutput, "MS%d, ", Nvtotal - nv); }
		fprintf(fMPDoutput, "\n");
		for (nf = 0; nf < Nftotal; nf++) {
			fprintf(fMPDoutput, "%d, ", nf);
			for (nv = 0; nv < Nvtotal; nv++) { 
				fprintf(fMPDoutput, "%lf, ", MPDoutput2nd[nv][nf][mrow]); }
			fprintf(fMPDoutput, "\n");		}
		fclose(fMPDoutput);		}
}
//																ここから印刷切れます↓
// 4.2.3  TotalCEFを、ファイバ長z毎に出力
void sub_routine_4_2_3(char dir[], int mrow, int noffV2F, int noffF2F, double ****TotalCEF) {
	FILE *fallCEFz;
	char str_sub[256];	str_sub[0] = '\0';	int nv;
	sprintf(str_sub, "%s\\%s_z_%dm.csv", dir, dirTotalCEF, mrow);
	fallCEFz = fopen(str_sub, "a+");	
	fprintf(fallCEFz, "%d, %d, ", noffV2F*OFFres_V2F, noffF2F*OFFres_F2F);
	for (nv = 0; nv < Nvtotal; nv++) {
		fprintf(fallCEFz, "%lf, ", TotalCEF[nv][noffV2F][noffF2F][mrow]);	}
	fprintf(fallCEFz, "\n");
	fclose(fallCEFz);			}

// 4.2.4  TotalCEFを、オフセット量 x1・x2 の組合せ毎に出力
void sub_routine_4_2_4(char dir[], int mrow, int noffV2F, int noffF2F, double ****TotalCEF) {
	FILE *fallCEFx1x2;
	char str_sub[256];	str_sub[0] = '\0';	int nv;
	sprintf(str_sub, "%s\\%s_x1=%dum_x2=%dum.csv", dir, dirTotalCEF, noffV2F*OFFres_V2F, noffF2F*OFFres_F2F);
	fallCEFx1x2 = fopen(str_sub, "a+");
	fprintf(fallCEFx1x2, "%d, %d, %d, ", noffV2F*OFFres_V2F, noffF2F*OFFres_F2F, mrow);
	for (nv = 0; nv < Nvtotal; nv++) {
		fprintf(fallCEFx1x2, "%lf, ", TotalCEF[nv][noffV2F][noffF2F][mrow]);
	}
	fprintf(fallCEFx1x2, "\n");
	fclose(fallCEFx1x2);
}

/*A.1. 第2種変形ベッセル関数 bessk (n, x) */
// 第1種変形Bessel関数（n=0）I0(x)
double bessi0(double x)
{
	double ax, ans, y;
	// Polynomial fit
	if ((ax = fabs(x)) < 3.75) {
		y = x / 3.75;
		y *= y;
		ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
			+ y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
	}
	else {
		y = 3.75 / ax;
		ans = (exp(ax) / sqrt(ax))*(0.39894228 + y * (0.1328592e-1
			+ y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
				+ y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
					+ y * 0.392377e-2))))))));
	}
	return ans;
}

// 第2種変形Bessel関数（n=0） K0(x)
double bessk0(double x)
{
	double bessi0(double x);
	double y, ans;
	//polynomial fit
	if (x <= 2) {
		y = x * x / 4.0;
		ans = (-log(x / 2.0)*bessi0(x)) + (-0.57721566 + y * (0.42278420
			+ y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2
				+ y * (0.10750e-3 + y * 0.74e-5))))));
	}
	else {
		y = 2.0 / x;
		ans = (exp(-x) / sqrt(x))*(1.25331414 + y * (-0.7832358e-1
			+ y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2
				+ y * (-0.251540e-2 + y * 0.53208e-3))))));
	}
	return ans;
}

// 第1種変形Bessel関数（n=1） I1(x)
double bessi1(double x)
{
	double ax, ans, y;
	if ((ax = fabs(x)) < 3.75) {
		y = x / 3.75;	   y *= y;
		ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
			+ y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
	}
	else {
		y = 3.75 / ax;
		ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1
			- y * 0.420059e-2));
		ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2
			+ y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
		ans *= (exp(ax) / sqrt(ax));
	}
	return x < 0.0 ? -ans : ans;
}

// 第2種変形Bessel関数（n=1） K1 (x)
double bessk1(double x) {
	double bessi1(double x);
	double y, ans;
	if (x <= 2.0) {
		y = x * x / 4.0;
		ans = (log(x / 2.0)*bessi1(x)) + (1.0 / x)*(1.0 + y * (0.15443144
			+ y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1
				+ y * (-0.110404e-2 + y * (-0.4686e-4)))))));
	}
	else {
		y = 2.0 / x;
		ans = (exp(-x) / sqrt(x))*(1.25331414 + y * (0.23498619
			+ y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2
				+ y * (0.325614e-2 + y * (-0.68245e-3)))))));
	}
	return ans;
}

// 第2種変形Bessel関数 Kn(x)
double bessk(int n, double x) {
	double bessk0(double x);
	double bessk1(double x);
	double bk, bkm, bkp, tox;   int j;
	if (n == 0) return bessk0(x);
	if (n == 1) return bessk1(x);
	if (n >= 2) {
		tox = 2.0 / x;
		bkm = bessk0(x);
		bk = bessk1(x);
		for (j = 1; j < n; j++) {
			bkp = bkm + j * tox*bk;
			bkm = bk;
			bk = bkp;
		}
		return bk;
	}
	else return 0;
}

/*A.2. 屈折率の波長微分関数dndl (n_lamda, mater) */
double dndl(double lamda, double n_lamda, int mater)
{
	if (mater == 0) // PMMA
		return ((0.02173e-4*lamda - 18.79107e-4)*n_lamda + (-0.03109e-4*lamda + 26.85035e-4))*1.0e9;
	else if (mater == 1) // CYTOP
		return ((0.02383e-5*lamda - 25.34898e-5)*n_lamda + (-0.02983e-5*lamda + 31.46375e-5))*1.0e9;
	else if (mater == 2) //Silica
		return ((0.02817e-5*lamda - 29.61218e-5)*n_lamda + (-0.03807e-5*lamda + 39.02595e-5))*1.0e9;
	else if (mater == 3) //PTCEMA
		return ((0.0076424e-4*lamda - 7.92635e-4)*n_lamda + (-0.0107534e-4*lamda + 11.10000e-4))*1.0e9;
	else printf(" Material number is not correct !\n"); exit(EXIT_FAILURE);
}



/*A.3. 係数行列要素計算関数 S_matrix (a, b, q, m, n, v, w, D)*/
/********************************************************************************************************************************************/
/// b[0]=S00, b[1]=S11, ...................... b[j]=Sjj ......  b[n]=Snn					// a[0]=___, a[1]=S01, a[2]=S12, ... a[j]=Sj-1,j ... a[n]=Sn-1,n					// a[0]=___, a[1]=S10, a[2]=S21, ... a[j]=Sj,j-1 ... a[n]=Sn,n-1
void S_matrix(double *a, double *b, double *q, int m, int n, double v, double w, double D) {
	if (m == 0) {
		b[0] = -(1.0 / 2.0) + (3.0*q[0] + 2.0*q[1])*(v*v / 60.0)*((D*D) / (n*n)) - (1.0 / 12.0)*(w*w)*((D*D) / (n*n));
		b[n] = -((2.0*n - 1.0) / 2.0) + ((5.0*n - 2.0)*q[n - 1] + 3.0*(5.0*n - 1.0)*q[n])*(v*v / 60.0)*((D*D) / (n*n)) - ((4.0*n - 1.0) / 12.0)*(w*w)*((D*D) / (n*n))
			- w * D*bessk(1, w*D) / bessk(0, w*D);
	}
	else {
		b[0] = 0.0;
		b[n] = -(1.0 - m * m)*((2.0*n - 1.0) / 2.0) + ((5.0*n - 2.0)*q[n - 1] + 3.0*(5.0*n - 1.0)*q[n])*(v*v / 60.0)*((D*D) / (n*n)) - ((4.0*n - 1.0) / 12.0)*(w*w)*((D*D) / (n*n))
			- (m*m)*((n - 1.0)*(n - 1.0))*log((double)n / ((double)n - 1.0)) - m * m - w * D*bessk(m - 1, w*D) / bessk(m, w*D) - m;
	}
	///----------------------------------------
	b[1] = -2.0*(1.0 - m * m) + (3.0*q[0] + 30 * q[1] + 7.0*q[2])*(v*v / 60.0)*((D*D) / (n*n)) - (2.0 / 3.0)*(w*w)*((D*D) / (n*n)) - (m*m)*4.0*log(2.0);			int j;
	///----------------------------------------
	for (j = 2; j < n; j++) {
		b[j] = -2.0*j*(1.0 - m * m) + ((5.0*j - 2.0)*q[j - 1] + 30.0*j*q[j] + (5.0*j + 2.0)*q[j + 1])*(v*v / 60.0)*((D*D) / (n*n)) - (2.0 / 3.0)*j*(w*w)*((D*D) / (n*n))
			- (m*m)*(((j - 1.0)*(j - 1.0))*log((double)j / ((double)j - 1.0)) + ((j + 1.0)*(j + 1.0))*log(((double)j + 1.0) / (double)j));
	}
	///----------------------------------------
	a[0] = 0.0;//未使用要素につき0を格納
	a[1] = (1.0 / 2.0) + (2 * q[0] + 3 * q[1])*(v*v / 60.0)*((D*D) / (n*n)) - (1.0 / 12.0)*(w*w)*((D*D) / (n*n)) - (m*m) / 2;
	///----------------------------------------
	for (j = 2; j <= n; j++) {
		a[j] = ((2.0*j - 1.0) / 2.0)*(1.0 - m * m) + ((5.0*j - 3.0)*q[j - 1] + (5.0*j - 2.0)*q[j])*(v*v / 60.0)*((D*D) / (n*n)) - ((2.0*j - 1.0) / 12.0)*(w*w)*((D*D) / (n*n))
			+ (m*m)*(j - 1.0)*j*log((double)j / ((double)j - 1.0));
	}
}
/********************************************************************************************************************************************/


/*A.4. 改訂コレスキー分解	mcholesky ( a, b, ML, MD, m, n )*/
void mcholesky(double *a, double *b, double *ML, double *MD, int m, int n) {
	int i;
	if (m == 0) {
		ML[0] = 0.0;//未使用要素につき0を格納
		MD[0] = b[0];
		for (i = 1; i <= n; i++) {
			MD[i] = b[i] - a[i] * a[i] / MD[i - 1];
			ML[i] = a[i] / MD[i - 1];
		}
	}
	else {
		ML[0] = 0.0;//未使用要素につき0を格納
		ML[1] = 0.0;
		MD[0] = 0.0;
		MD[1] = b[1];
		for (i = 2; i <= n; i++) {
			MD[i] = b[i] - a[i] * a[i] / MD[i - 1];
			ML[i] = a[i] / MD[i - 1];
		}
	}
}


/*A.5. 改訂コレスキー分解法による求解	mcholesky_sol ( a, b, ML, MD, m, n )*/
void mcholesky_sol(double *ML, double *MD, double *R, int m, int n) {
	int i;
	//「Ly=R0」を解く
	if (m == 0) {
		for (i = 1; i <= n; i++) {
			R[i] = R[i] - R[i - 1] * ML[i];
		}
		//「(D(LT))R1=y」を「(LT)R1=(D-1)y=y'」に変える
		for (i = 1; i <= n; i++) {
			R[i] = R[i] / MD[i];
		}
		//「 (LT)R1=y'」を解く
		for (i = n - 1; i >= 0; i--) {
			R[i] = R[i] - ML[i + 1] * R[i + 1];
		}
	}

	else {
		for (i = 2; i <= n; i++) {
			R[i] = R[i] - R[i - 1] * ML[i];
		}
		//「(D(LT))R1=y」を「(LT)R1=(D-1)y=y'」に変える
		for (i = 2; i <= n; i++) {
			R[i] = R[i] / MD[i];
		}
		//「 (LT)R1=y'」を解く
		for (i = n - 1; i >= 1; i--) {
			R[i] = R[i] - ML[i + 1] * R[i + 1];
		}
	}
}

/*A.6. 逆べき乗法の初期ベクトル計算	R0 ( MD, R, m, n )*/
/* 対角行列Dの成分が最大となる要素だけ1であるようなベクトルを選定*/
void R0(double *MD, double *R, int m, int n) {
	int i, j = 1;
	if (m == 0) {
		for (i = 0; i <= n - 1; i++) {
			if (fabs(MD[i + 1]) < fabs(MD[j])) { j = i + 1; }
		}
	}
	else {
		for (i = 1; i <= n - 1; i++) {
			if (fabs(MD[i + 1]) < fabs(MD[j])) { j = i + 1; }
		}
	}
	//R0の初期値の代入
	for (i = 0; i <= n; i++) {
		if (i == j) { R[i] = 1.0; }
		else { R[i] = 0.0; }
	}
}

/*A.7. 逆べき乗法の解ベクトル規格化	R_norm ( R, n )*/
void R_norm(double *R, int n) {
	int i;	double s = 0;
	for (i = 0; i <= n; i++) { s = s + R[i] * R[i]; }			//行列要素の２乗和
	if (s != 0) {
		for (i = 0; i <= n; i++) { R[i] = R[i] / sqrt(s); }
	}
}

/*A.8. 固有値計算（Rayleigh quotient）Eigen ( R, a, b, m, n )*/// Rベクトルを規格化しているため内積は1
double Eigen(double *R, double *a,
	double *b, int n, int m) {
	int i;	double s = 0;
	if (m == 0) {
		s = (R[0] * b[0] + R[1] * a[1])*R[0];
		for (i = 1; i < n; i++) {
			s += (R[i - 1] * a[i] + R[i] * b[i] + R[i + 1] * a[i + 1])*R[i];
		}
		s += (R[n - 1] * a[n] + R[n] * b[n])*R[n];
		return s;
	}
	else {
		s = (R[1] * b[1] + R[2] * a[2])*R[1];
		for (i = 2; i < n; i++) {
			s += (R[i - 1] * a[i] + R[i] * b[i] + R[i + 1] * a[i + 1])*R[i];
		}
		s += (R[n - 1] * a[n] + R[n] * b[n])*R[n];
		return s;
	}
}

/*A.9. 群遅延計算用関数	dbdk_bunbo ( R, D, w, m, n )*/
/* 入力パラメータ（横方向電場分布，コア径，伝搬定数，要素分割数）に対するm次モード群遅延計算式の分母分子を計算する */

// 横方向電場成分 R[0]~R[n], 規格化伝搬定数 w, 規格化コア径 D, 方位角モード次数 m, 分割数 n
double dbdk_bunbo(double *R, double D, double w, int m, int n) {
	int i;	double s = 0;
	for (i = 0; i <= n - 1; i++) {
		s = s + (1.0 / 12.0) * ((D*D) / (n*n)) * ((double)(4 * i + 1)*R[i] * R[i] + 2.0*(double)(2 * i + 1)*R[i] * R[i + 1] + (double)(4 * i + 3)*R[i + 1] * R[i + 1]);
	}
	if (m == 0) {
		return s + ((bessk(1, w*D)*bessk(1, w*D) / (bessk(0, w*D)*bessk(0, w*D))) - 1.0) * ((D*D)*(R[n] * R[n]) / 2.0);
	}
	else {
		return s + ((bessk(m - 1, w*D)*bessk(m + 1, w*D) / (bessk(m, w*D)*bessk(m, w*D))) - 1.0) * ((D*D)*(R[n] * R[n]) / 2.0);
	}
}

/*A.10. 群遅延計算用関数	dbdk_bunshi ( R,qg, D, w, m, n )*/
/* 入力パラメータ（横方向電場分布，コア径，伝搬定数，要素分割数）に対するm次モード群遅延計算式の分母分子を計算する */

// 横方向電場成分 R[0]~R[n], 規格化伝搬定数 w, 規格化コア径 D, 方位角モード次数 m, 分割数 n
// 屈折率分散パラメータ qg[0]~qg[n] ( = n*(d(kn)/dk) )
double dbdk_bunshi(double *R, double *qg, double D, double w, int m, int n)
{
	int i;		double s = 0;
	for (i = 0; i <= n - 1; i++) {
		s = s + (1.0 / 12.0) * ((D*D) / (n*n))* (((double)(3 * i) + 3.0 / 5.0)*qg[i] * R[i] * R[i] + ((double)i + 2.0 / 5.0)*(2.0*qg[i] * R[i + 1] + qg[i + 1] * R[i])*R[i] + ((double)i + 3.0 / 5.0)*(qg[i] * R[i + 1] + 2.0*qg[i + 1] * R[i])*R[i + 1] + ((double)(3 * i) + 12.0 / 5.0)*qg[i + 1] * R[i + 1] * R[i + 1]);
	}
	if (m == 0) {
		return s + qg[n] * ((bessk(1, w*D)*bessk(1, w*D) / (bessk(0, w*D)*bessk(0, w*D))) - 1.0) * ((D*D)*(R[n] * R[n]) / 2.0);
	}
	else {
		return s + qg[n] * ((bessk(m - 1, w*D)*bessk(m + 1, w*D) / (bessk(m, w*D)*bessk(m, w*D))) - 1.0) * ((D*D)*(R[n] * R[n]) / 2.0);
	}
}

/* A.11. 畳み込み積分関数 convolution ( xmax, P, LT, pulse) */
double *convolution(int xmax, double *P, int LT, double *pulse) {
	int i, j, xL;
	double *R;
	xL = sizeof(double)*(xmax + LT + 1);
	R = (double*)malloc(xL);
	printf("xmax = %d\n" "LT = %d\n" "xL = %d\n", xmax, LT, xL);
	//for ( i = 0 ; i <= xmax ; i++ ){ printf ( "%f\n", P[i] );}
	//for ( i = 0 ; i <= LT ; i++ ){ printf ( "%f\n", pulse[i] );}
	for (j = 0; j <= xmax + LT; j++) { R[j] = 0; }
	for (i = 0; i < xmax; i++) {
		for (j = 0; j <= LT; j++) {
			R[i + j] = R[i + j] + P[i] * pulse[j];
			//printf ( "%f\n", R[i+j] );
		}
	}				//for ( i = 0; i <= xmax+LT+1 ; i++ ){ printf ( "%f\n", R[i] );}
	return R;
}

/* A.12. 整数ベクトル領域確保用関数 dvector (i, j) */
int *dintvector(int i, int j) {
	int *a;
	if ((a = (int *)malloc((j - i + 1) * sizeof(int))) == NULL) {
		printf("Memory cannot be allocated !\n"); exit(1);
	}
	return (a - i);
}

/* A.13. 整数ベクトル領域解放用関数 free_dvector (a, i) */
void free_dintvector(int *a, int i) {
	free((void*)(a + i));
}


/* A.14. 整数ベクトル領域確保用関数 dvector (i, j) */
double *drealvector(int i, int j) {
	double *a;
	if ((a = (double *)malloc((j - i + 1) * sizeof(double))) == NULL) {
		printf("Memory cannot be allocated !\n"); exit(1);
	}
	return (a - i);
}

/* A.15. 整数ベクトル領域解放用関数 free_dvector (a, i) */
void free_drealvector(double *a, int i) {
	free((void*)(a + i));
}


/* A.16. 実数行列領域確保用関数 dmatrix ( nr1, nr2, nl1, nl2 ) */
double **dmatrix(int nr1, int nr2, int nl1, int nl2) {
	// nrow: 行の数, ncol: 列の数
	double **a;
	int i, nrow, ncol;
	nrow = nr2 - nr1 + 1;
	ncol = nl2 - nl1 + 1;
	/* 行の確保 */
	if ((a = (double **)malloc(nrow * sizeof(double*))) == NULL) {
		printf("Memory cannot be allocated !\n"); exit(1);
	}
	a = a - nr1; // 行をずらす
				 /* 列の確保 */
	for (i = nr1; i <= nr2; i++) a[i] = (double *)malloc(ncol * sizeof(double));
	for (i = nr1; i <= nr2; i++) a[i] = a[i] - nl1;	// 列をずらす
	return (a);
}

/* A.17.実数行列領域解放用関数 free_dmatrix ( a, nr1, nr2, nl1, nl2 ) */
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2) {
	int i;
	for (i = nr1; i <= nr2; i++) free((void*)(a[i] + nl1));
	free((void*)(a + nr1));
}


/* A.18. 整数ベクトル初期化関数 init_vector ( a, nr1, nr2 ) */
void init_intvector(int *a, int nr1, int nr2) {
	int i;
	for (i = nr1; i <= nr2; i++) { a[i] = 0; }
}


/* A.19. 整数ベクトル初期化関数 init_vector ( a, nr1, nr2 ) */
void init_realvector(double *a, int nr1, int nr2) {
	int i;
	for (i = nr1; i <= nr2; i++) { a[i] = 0.0; }
}

/* A.21. エラー関数 -------(A.22とA.23に用いる)-----------*/
void nrerror(char error_text[]) {
	fprintf(stderr, "Numerical Recipes run=time error///\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}


/* A.22. float型の2次元配列 matrix(long nrl, long nrh, long ncl, long nch) */
float **matrix(long nrl, long nrh, long ncl, long nch) {
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float **m;
	m = (float **)malloc((size_t)((nrow + NR_END) * sizeof(float*)));
	// if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;		  m -= nrl;
	m[nrl] = (float *)malloc((size_t)((nrow*ncol + NR_END) * sizeof(float*)));
	// if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;	  m[nrl] -= ncl;
	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
	return m;
}

/* A.23. float型2次元配列の解放 matrix(long nrl, long nrh, long ncl, long nch) */
void free_matrix(float **a, long nrl, long nrh, long ncl, long nch) {
	free((FREE_ARG)(a[nrl] + ncl - NR_END));
	free((FREE_ARG)(a + nrl - NR_END));
}

/* A.24. float型2次元配列の解放 f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) */
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
	long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	double ***t;
	/* allocate pointers to pointers to rows */
	t = (double ***)malloc((size_t)((nrow + NR_END) * sizeof(double**)));
	//if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;
	/* allocate pointers to rows and set pointers to them */
	t[nrl] = (double **)malloc((size_t)((nrow*ncol + NR_END) * sizeof(double*)));
	//if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	/* allocate rows and set pointers to them */
	t[nrl][ncl] = (double *)malloc((size_t)((nrow*ncol*ndep + NR_END) * sizeof(double)));
	//if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	for (j = ncl + 1;j <= nch;j++) t[nrl][j] = t[nrl][j - 1] + ndep;
	for (i = nrl + 1;i <= nrh;i++) {
		t[i] = t[i - 1] + ncol;
		t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
		for (j = ncl + 1;j <= nch;j++) t[i][j] = t[i][j - 1] + ndep;
	}
	/* return pointer to array of pointers to rows */
	return t;
}

/* A.24. float型2次元配列の解放 f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) */
void free_d3tensor(double ***a, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
	/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG)(a[nrl][ncl] + ndl - NR_END));
	free((FREE_ARG)(a[nrl] + ncl - NR_END));
	free((FREE_ARG)(a + nrl - NR_END));
}

/* A.25. 実数行列初期化関数 init_dmatrix ( a, nr1, nr2, nl1, nl2 ) */
void init_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2) {
	for (int i = nr1; i <= nr2; i++) {
		for (int j = nl1; j <= nl2; j++) { a[i][j] = 0.0; }
	}
}

/* A.26. 実数行列初期化関数 init_d3tensor ( a, nr1, nr2, nl1, nl2, np1, np2 ) */
void init_d3tensor(double ***a, int nr1, int nr2, int nl1, int nl2, int np1, int np2) {
	for (int i1 = nr1; i1 <= nr2; i1++) {
		for (int i2 = nl1; i2 <= nl2; i2++) { 
			for(int i3 = np1; i3 <= np2; i3++) { a[i1][i2][i3] = 0.0; }
}}}

/* A.27. 実数行列初期化関数 init_matrix ( a, nr1, nr2, nl1, nl2 ) */
void init_matrix(float **a, int nr1, int nr2, int nl1, int nl2) {
	for (int i = nr1; i <= nr2; i++) {
		for (int j = nl1; j <= nl2; j++) { a[i][j] = 0.0; }
	}
}

/* A.28. 4次元配列 ( a, nr1, nr2, nl1, nl2 ) */
double ****d4tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long ntl, long nth) {
	long j1, j2, j3; long nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1, ntime = nth - ntl + 1;
	double ****t;
	/* allocate pointers to pointers to rows */
	t = (double ****)malloc((size_t)((nrow + NR_END) * sizeof(double***)));
	//if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;
	/* allocate pointers to rows and set pointers to them */
	t[nrl] = (double ***)malloc((size_t)((nrow*ncol + NR_END) * sizeof(double**)));
	//if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	/* allocate rows and set pointers to them */
	t[nrl][ncl] = (double **)malloc((size_t)((nrow*ncol*ndep + NR_END) * sizeof(double*)));
	//if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	/* allocate rows and set pointers to them */
	t[nrl][ncl][ndl] = (double *)malloc((size_t)((nrow*ncol*ndep*ntime + NR_END) * sizeof(double)));
	t[nrl][ncl][ndl] += NR_END;
	t[nrl][ncl][ndl] -= ntl;
	for (j3 = ndl + 1; j3 <= ndh; j3++) t[nrl][ncl][j3] = t[nrl][ncl][j3 - 1] + ntime;
	for (j2 = ncl + 1; j2 <= nch; j2++) {
		t[nrl][j2] = t[nrl][j2 - 1] + ndep;
		t[nrl][j2][ndl] = t[nrl][j2 - 1][ndl] + ndep * ntime;
		for (j3 = ndl + 1; j3 <= ndh; j3++) t[nrl][j2][j3] = t[nrl][j2][j3 - 1] + ntime;	
	}
	for (j1 = nrl + 1; j1 <= nrh; j1++) {
		t[j1] = t[j1 - 1] + ncol;
		t[j1][ncl] = t[j1 - 1][ncl] + ncol * ndep;
		t[j1][ncl][ndl] = t[j1 - 1][ncl][ndl] + ncol * ndep * ntime;
		for (j3 = ndl + 1; j3 <= ndh; j3++) t[j1][ncl][j3] = t[j1][ncl][j3 - 1] + ntime;
		for (j2 = ndl + 1; j2 <= ndh; j2++) {
			t[j1][j2] = t[j1][j2 - 1] + ndep;
			t[j1][j2][ndl] = t[j1][j2 - 1][ndl] + ndep * ntime;
			for (j3 = ndl + 1; j3 <= ndh; j3++) t[j1][j2][j3] = t[j1][j2][j3 - 1] + ntime;
		}
	}
	return t;
}
void free_d4tensor(double ****a, long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long ntl, long nth) {
	free((FREE_ARG)(a[nrl][ncl][ndl] + ntl - NR_END));
	free((FREE_ARG)(a[nrl][ncl] + ndl - NR_END));
	free((FREE_ARG)(a[nrl] + ncl - NR_END));
	free((FREE_ARG)(a + nrl - NR_END));
}
void init_d4tensor(double ****a, long nr1, long nr2, long nc1, long nc2, long nd1, long nd2, long nt1, long nt2) {
	for (int i1 = nr1; i1 <= nr2; i1++) {
		for (int i2 = nc1; i2 <= nc2; i2++) {
			for (int i3 = nd1; i3 <= nd2; i3++) {
				for (int i4 = nt1; i4 <= nt2; i4++) { a[i1][i2][i3][i4] = 0.0; }}}}
}

/* d3tensor の第2、第3成分の初期化*/
void init_d3tensor_nfnz(double ***a, int nr, int nl1, int nl2, int np1, int np2){
	for (int i2 = nl1; i2 <= nl2; i2++) { 
		for(int i3 = np1; i3 <= np2; i3++) { a[nr][i2][i3] = 0.0; }
}}

//																ここから印刷切れます↓
void mkdir(char dirname[]) {
	struct stat statBuf;
	if (stat(dirname, &statBuf) != 0){
		if (_mkdir(dirname) != 0) {
			printf("ディレクトリ %s の作成に失敗しました。\n", dirname);
			system("pause");	}}}