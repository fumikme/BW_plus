/** 更新履歴 **/
/********* 2019/04/09  **********///											 ↓ここから印刷切れます
/* FEM_LPmode_VCSELinputを、コア半径10umのステップ型分布ファイバの分布として計算したFEM_resultから算出している点に注意 */

/*******************************************************************************************************/
/*******************************                              ******************************************/
/*******************************         関心毎の分離         ******************************************/
/*******************************                              ******************************************/
/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
//                                                                                                     //
/*-*-                     モジュール化してそれぞれ10行くらいにしておくと，                        -*-*-*/
/*-*-             実装する時にも分離した小さい単位のロジックに思考を集中できるので、              -*-*-*/
/*-*-                              コードが書きやすくなります。                                   -*-*-*/
//                                                                                                     //
/*******************************************************************************************************/

#define _USE_MATH_DEFINES
#define NR_END 1
#define FREE_ARG char*
#pragma warning(disable:4996) //fopen関数を許可
#pragma warning(disable:0266)

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;


/* 諸定数の定義 */
#define PI 3.141592653 // 円周率
#define C 2.99792458e8 // 光速
#define Epsilon0 8.8541878e-12 // 真空中の誘電率
#define Mu0 PI*4.0e-7 // 真空中の透磁率

#define num_mode_vin 6

/* 関数副プログラムの読み込み *//* 関数副プログラムの読み込み */
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
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void init_matrix(float **a, int nr1, int nr2, int nl1, int nl2);
// 2次元配列（double型）
double  **dmatrix(int nr1, int nr2, int nl1, int nl2);
void free_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
void init_dmatrix(double **a, int nr1, int nr2, int nl1, int nl2);
// 3次元配列
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void init_d3tensor(double ***a, int nr1, int nr2, int nl1, int nl2, int np1, int np2);

/* エラー出力関数の読み込み */
void nrerror(char error_text[]);

/* 0. パラメータの定義 */
// m: モード次数( TE&TM ~ 1, EH ~ n+1, HE ~ n-1), l: 動径方向モード次数, n: 方位角モード次数
// mater: 材料ID, profile: 屈折率分布入力方法ID (0: べき乗分布, 1: 測定分布入力)
// NLP: LPモード数, NLP0: 0次方位角モードLPモード数, Ntotal: 全モード数
// N: 動径座標rのコア内空間ステップ数, Nclad: 解析領域空間ステップ数, Nbeta: 伝搬定数の分割数
// Nwkb: WKB近似により求めたLPモード数, Ntmin: 最小群遅延時間ステップ数, Ntmax: 最大群遅延時間ステップ数
// lamda: 波長, k: 波数, omega: 角周波数, delta: 比屈折率, NA: 開口数, aa: 動径座標規格化サイズ
// A: コア半径, AA: 解析領域半径, g: べき乗分布指数, n1: コア中心の屈折率, n0: クラッドの屈折率
// L: ファイバ長, Tv: インパルス応答時間軸刻み幅, dr: 動径座標刻み幅,
// v: 規格化周波数, w: 規格化伝搬定数, D: 規格化コア径
// tau: 群速度, beta: 伝搬定数β, dbeta: 伝搬定数刻み幅
// a: 係数行列Sの副対角要素格納配列, b: 係数行列Sの対角要素格納配列
// GI: コア内屈折率, q: 規格化コア内屈折率, qg: コア内屈折率分散パラメータ
// R: 規格化横方向電場成分, Rb: 逆べき乗法用入れ子配列, Rinf: 解析領域
// ML: LDU分解係数行列のL行列副対角要素, MD: LDU分解係数行列のD行列対角要素
// de, df: 逆べき乗法における収束判別パラメータ, eig: 固有値
// taumin: 最小単位群遅延, taumax: 最大単位群遅延
// fmin: 評価最小周波数, fmax: 評価最大周波数, dfrq: 評価周波数分解能
// Hw: 伝達関数要素入れ子変数 (ReHw: 実部, ImHw: 虚部), bw: -3dB帯域幅
// GI: 屈折率分布, q: 屈折率波長微分, qg: n*d(kn)/dk
// Rlp: 全LPモードのモードフィールド, Mbeta: 各モードの伝搬定数
// r0: オフセット（励振位置），w0: ガウシアンビーム半径，dx: 重なり積分の空間分解能, Nxy: 解析領域の分割数


int main(void)
{
	FILE   *fp, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8, *fr, *fpulse, *fopulse, *fmpd;
	char   fppath[128], fp2path[128], fp3path[128], fp4path[128], fp5path[128], fp6path[128], fp7path[128];
	char   fp8path[128], frpath[128], fpulsepath[128], fopulsepath[128], fmpdpath[128];
	          fppath[0]   =	 fp2path[0] = fp3path[0] = fp4path[0] = fp5path[0] = fp6path[0] = fp7path[0] =
	          fp8path[0]  =  frpath[0] = fpulsepath[0] = fopulsepath[0] = fmpdpath[0] = '\0' ;
	int    myu, x, y, i, j, jmax, count;
	int    m, l, NLP, NLP0, Ntotal, N, Nclad, Nbeta, mater, profile, Nwkb, Ntmin, Ntmax, Nf, launch;
	int    Nxy, xmax, LT, xL;
	double lamda, k, omega, A, AA, g, n0, n1, dr, L, Tv;
	double r0, w0, dx, dy;
	double delta, NA, aa, v, w, D;
	double tau, beta, dbeta, bb, eps1, eps2, sum, Rinf, de, df, eig;
	double taumin, taumax, fmin, fmax, dfrq, Hw, ReHw, ImHw, bw;					//	int    nr;
	double yy, Rxy, Rinxy;
	double Em_evin, Em_odin, Em_ev, Em_od, Am_evev, Am_evod, Am_odev, Am_odod;
	double ncav, noxi, Avin;
	double *GI, *pulse, *opulse, *q, *qg, *R, *R2, *Rb, *a, *b, *ML, *MD, *Mtau, *P, *M, *Mbeta, *MPD;
	double *Mbetain;
	double **Rlp, **MPD2d, **Rinlp;;
	double tauin, betain, Rinfin, eigin, betain_devided_by_k, data;
	double win, xx1, xx2, rr1, rr2;
	int    *modem, *modemin, *model, *modelin;
	int    Nvin, couple, min, lin, nrr1, nrr2, max;
	int    OFFres, OFFrange;

	char   trash[65536] = "\0";
	char   directory[128] = "FILE_for_IO";

	/** 1. 入力ファイルの読み込み **/
	if ((fp = fopen("FEM_input_kob.csv", "r")) != NULL) {
		char s1[256], s2[256];
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &A);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &AA);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &profile);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &g);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &n1);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &n0);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &launch);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &couple);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &r0);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &OFFrange);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &OFFres);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &Avin);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &w0);			// ガウス分布（シングルモードレーザ）を用いて励振
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &dr);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &dx);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &mater);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &lamda);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &L);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &Tv);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &fmin);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &fmax);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &Nf);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &aa);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &Nbeta);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &eps1);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &eps2);
		fscanf(fp, "%[^,], %[^,], %d\n", s1, s2, &jmax);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &noxi);
		fscanf(fp, "%[^,], %[^,], %lf\n", s1, s2, &ncav);
		// fscanf(fp, "%[^,]\n", s1);
	}
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }

	N = (int)(1000.0*A / dr); Nclad = (int)(1000.0*(AA - A) / dr);				//-----------Ncladに注意
	Nvin = (int)(1000.0*Avin / dr);


	fmin = fmin * 1.0e-3; fmax = fmax * 1.0e-3; dfrq = (fmax - fmin) / (double)Nf;

	/************************   ここが非常に重要   ****************************/
	Nxy = (int)(3 * 1000 * A / dx); //限定モード励振時の解析領域に対応する分割数 Nxy*dx=解析領域			
	/**************************************************************************/

	lamda = lamda * 1.0e-9;  A = A * 1.0e-6; dr = dr * 1.0e-9; AA = AA * 1.0e-6;
	Avin = Avin * 1.0e-6;
	r0 = r0 * 1.0e-6; w0 = w0 * 1.0e-6, dx = dx * 1.0e-9; dy = dx;	// 単位変換（THz(1/ps)）			
	fclose(fp);

	/* 測定プロファイル読み込み */
	GI = drealvector(0, N); init_realvector(GI, 0, N);
	sprintf(frpath, "%s/profile.csv", directory);
	if (profile == 1) {
		if ((fr = fopen(frpath, "r")) != NULL) {
			for (j = 0; j <= N; j++) { fscanf(fr, "%lf,", &GI[j]); }
		}
		else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }
		n1 = GI[0]; n0 = GI[N]; fclose(fr);
	}

	/* 入射波形読み込み */
	pulse = drealvector(0, 1999); init_realvector(pulse, 0, 1999); //入射波形は2000ステップで入力
	sprintf(fpulsepath, "%s/input pulse.csv", directory);
	if ((fpulse = fopen(fpulsepath, "r")) != NULL) {
		for (j = 0; j <= 1999; j++) { fscanf(fpulse, "%lf,", &pulse[j]); }
	}
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }
	fclose(fpulse);


	/** 2. 各種定数の設定 **/
	k = 2.0*PI / lamda; omega = 2.0*PI*C / lamda;
	delta = (n1*n1 - n0 * n0) / (2.0*n1*n1); NA = sqrt(n1*n1 - n0 * n0);
	v = k * aa*n1*sqrt(2.0*delta); D = A / aa;

	Nwkb = int((1.0 / 4.0)*(g / (g + 2.0))*(k*k)*(n1*n1)*delta*(A*A));
	bb = -1; dbeta = k * (n1 - n0) / (double)Nbeta;//初期化
	NLP = NLP0 = Ntotal = 0;//初期化

							/* 配列の記憶領域確保および初期化 */
	q = drealvector(0, N); init_realvector(q, 0, N);
	qg = drealvector(0, N); init_realvector(qg, 0, N);
	R = drealvector(0, N); init_realvector(R, 0, N);
	R2 = drealvector(0, N); init_realvector(R2, 0, N);
	Rb = drealvector(0, N); init_realvector(Rb, 0, N);
	a = drealvector(0, N); init_realvector(a, 0, N);
	b = drealvector(0, N); init_realvector(b, 0, N);
	ML = drealvector(0, N); init_realvector(ML, 0, N);
	MD = drealvector(0, N); init_realvector(MD, 0, N);
	M = drealvector(0, Nf); init_realvector(M, 0, Nf);
	Mtau = drealvector(0, 10 * Nwkb); init_realvector(Mtau, 0, 10 * Nwkb);
	modem = dintvector(0, 10 * Nwkb); init_intvector(modem, 0, 10 * Nwkb);
	model = dintvector(0, 10 * Nwkb); init_intvector(model, 0, 10 * Nwkb);
	modelin = dintvector(0, 10 * Nwkb); init_intvector(modelin, 0, 10 * Nwkb);
	Mbeta = drealvector(0, 10 * Nwkb); init_realvector(Mbeta, 0, 10 * Nwkb);
	MPD = drealvector(0, 10 * Nwkb); init_realvector(MPD, 0, 10 * Nwkb);

	#define sizeofRLP 500
	Rlp = dmatrix(0, sizeofRLP, 0, sizeofRLP);
	for (int i = 0; i <= sizeofRLP; i++) {
		for (int j = 0; j <= sizeofRLP; j++) { Rlp[i][j] = 0.0; }
	}

	/* 屈折率分布の設定 */
	/* べき乗プロファイル*/
	if (profile == 0) {
		for (j = 0; j <= N; j++) { GI[j] = n1 * sqrt(1.0 - 2.0*delta*pow(((double)j / (double)N), g)); }
	}
	for (j = 0; j <= N; j++) { q[j] = (GI[j] * GI[j] - n0 * n0) / (n1*n1 - n0 * n0); }
	for (j = 0; j <= N; j++) { qg[j] = GI[j] * GI[j] - (lamda*GI[j] * dndl(lamda*1.0e9, GI[j], mater)) / (1 - (lamda / GI[j])*dndl(lamda*1.0e9, GI[j], mater)); }


	/** 3. 評価条件の出力 **/
	sprintf(fp2path, "%s/FEM_setting.csv", directory);
	if ((fp2 = fopen(fp2path, "w")) != NULL) {
		fprintf(fp2, "Material, mater, %d\n", mater);
		fprintf(fp2, "Wavelength, λ, %lf,nm\n", lamda*1e9);
		fprintf(fp2, "Fiber length, L, %lf,m\n", L);
		fprintf(fp2, "Core radius, A, %lf,μm\n", A*1e6);
		fprintf(fp2, "Analysis region in radial axis, AA,%lf,μm\n", AA*1e6);
		fprintf(fp2, "Index exponent, g, %lf\n", g);
		fprintf(fp2, "Refractive index at the core center, n1,%lf\n", n1);
		fprintf(fp2, "Refractive index in the cladding, n0,%lf\n", n0);
		fprintf(fp2, "Relative refractive index, Δ,%lf\n", delta);
		fprintf(fp2, "Numerical aperture, NA, %lf\n", NA);
		fprintf(fp2, "Step size of the elements, dr, %lf, nm\n", dr*1e9);
		fprintf(fp2, "Step size of the pulse waveform, Tv, %lf, ps\n", Tv);
		fprintf(fp2, "Minimum evaluated frequency,fminx,%e,GHz\n", fmin*1.0e3);
		fprintf(fp2, "Maximum evaluated frequency,fmax,%e,GHz\n", fmax*1.0e3);
		fprintf(fp2, "Step size of evaluated frequency,dfrq,%e,GHz\n", dfrq*1.0e3);
		fprintf(fp2, "Step size of propagation constant, dβ,%lf\n", dbeta);
		fprintf(fp2, "Partition number of propagation constant, Nβ,%d\n", Nbeta);
		fprintf(fp2, "Partition number of fiber core radius, N,%d\n", N);
		fprintf(fp2, "Partition number of fiber cladding, Nclad,%d\n", Nclad);
		fprintf(fp2, "Maximum allowable error for convergence solution vector, eps1,%lf\n", eps1);
		fprintf(fp2, "Maximum allowable error of zero eigen value, eps2,%lf\n", eps2);
		fprintf(fp2, "Maximum number of iterations in inverse power method, jmax,%d\n", jmax);
		printf("Material no.: %d\n", mater);
		printf("Wavelength λ: %lfnm\n", lamda*1.0e9);
		printf("Fiber length L: %lfm\n", L);
		printf("Core radius A: %lfμm\n", A*1.0e6);
		printf("Analysis region AA: %lf μm\n", AA*1e6);
		printf("Index exponent g: %lf\n", g);
		printf("Refractive index at the core center n1: %lf\n", n1);
		printf("Refractive index in the cladding n0: %lf\n", n0);
		printf("Relative refractive index Δ: %lf\n", delta);
		printf("Numerical aperture NA: %lf\n", NA);
		printf("Step size of the elements dr: %lfnm\n", dr*1.0e9);
		printf("Partition number of fiber core radius N: %d\n", N);
		printf("Partition number of fiber cladding Nclad: %d\n", Nclad);
		printf("Step size of propagation constants dβ: %lf\n", dbeta);
		printf("Partition number of propagation constants Nβ: %d\n", Nbeta);
		printf("Maximum allowable error for convergence solution vector eps1: %lf\n", eps1);
		printf("Maximum allowable error of zero eigen value eps2: %lf\n", eps2);
		printf("Maximum number of iterations in inverse power method, jmax,%d\n", jmax);
	}
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }

	sprintf(fp3path, "%s/FEM_profile.csv", directory);
	if ((fp3 = fopen(fp3path, "w")) != NULL) {
		fprintf(fp3, "r [μm], n, q, qg\n");
		for (j = 0; j <= N; j++) {
			fprintf(fp3, "%lf, %lf, %lf, %lf\n", A*1.0e6* (double)j / (double)N, GI[j], q[j], qg[j]);
		}
	}
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }
	fclose(fp3);

	if ((fp4 = fopen("FEM_result.csv", "w")) != NULL) {
		printf("m, l, tau, beta, ne, eig, R[N]\n");
		fprintf(fp4, "m,l,tau,beta,ne,R_infinite,eig,r=0_um,");
		for (j = 1; j <= N; j++) { fprintf(fp4, "%lf,", (double)j*dr*1.0e6); } fprintf(fp4, "\n");
	}
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }

	sprintf(fp5path, "%s/FEM_impulse_response.csv", directory);
	if ((fp5 = fopen(fp5path, "w")) != NULL) { ; }
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }

	sprintf(fp6path, "%s/FEM_frequency-respnse.csv", directory);
	if ((fp6 = fopen(fp6path, "w")) != NULL) {
		for (y = 0; y <= Nf; y++) { fprintf(fp6, "%lf,", (fmin*1.0e3) + (double)y*(dfrq*1.0e3)); } fprintf(fp6, "\n");
	}
	else { printf("The file cannot be opened !\n"); exit(EXIT_FAILURE); }

	sprintf(fopulsepath, "%s/FEM_output pulse.csv", directory);
	if ((fopulse = fopen(fopulsepath, "w")) != NULL) { ; }
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }

	sprintf(fmpdpath, "%s/FEM_mode power distribution.csv", directory);
	if ((fmpd = fopen(fmpdpath, "w")) != NULL) { ; }
	else { printf(" U cannot open the file !\n"); exit(EXIT_FAILURE); }
	

	/** 4. 固有ベクトルおよび固有値の計算 **/
	for (m = 0; ; m++) {
		l = 1; beta = k * n1;

		/****************************************************************************************/
		for (; ; ) {
			/* 係数行列計算および改訂コレスキー分解 */
			w = aa * sqrt((beta*beta) - (k*k)*(n0*n0));
			S_matrix(a, b, q, m, N, v, w, D);
			mcholesky(a, b, ML, MD, m, N);
			/* 初期ベクトル R0 の付与 */
			R0(MD, R, m, N);
			/*  連立一次方程式 SR=(LDL)R=bR の反復評価 */
			for (j = 0; ; j++) {
				for (i = 0; i <= N; i++) { Rb[i] = R[i]; }
				mcholesky_sol(ML, MD, R, m, N);
				R_norm(R, N);
				/* 収束判定 */
				de = 0, df = 0;
				for (i = 0; i <= N; i++) {
					de = de + (Rb[i] - R[i])*(Rb[i] - R[i]);
					df = df + (Rb[i] + R[i])*(Rb[i] + R[i]);
				}

				if (de < eps1 || df < eps1) break;
				if (j >= jmax) goto nextin;
				// ① RとRbの成分差deが0に漸近すれば収束（break）．
				// ② RとRbの成分和dfが0に漸近すれば中止（break）．
				// ③ 反復回数が上限値 jmax を超えたらβを変更して再計算（go to next）．
			}

			/* 固有値の計算 */
			eig = Eigen(R, a, b, N, m);

			/* 固有値の妥当性評価 */
			// 「収束固有値 eig が前回値 bb と同値」であればβを変えて初めから再計算
			if (eig == bb) {
				dbeta = k * (n1 - n0) / (double)Nbeta;
				beta = beta - 1.0*dbeta;
				continue;
			}

			//  ①「0< 収束固有値 eig < 0.00001」であれば零固有値として採用
			if (0.0 < eig && eig < eps2) {
				/* 横方向電場成分Rの規格化 （パワーを1Wとする）*/
				sum = 0.0;
				for (j = 0; j < N; j++) { sum = sum + R[j] * R[j] * (j*dr)*dr; }
				for (j = 0; j < Nclad; j++) { sum = sum + R[N] * (bessk(m, w*(j*dr + A)) / bessk(m, w*A))*R[N] * (bessk(m, w*(j*dr + A)) / bessk(m, w*A))*(j*dr + A)*dr; }
				for (j = 0; j <= N; j++) { R2[j] = R[j] * sqrt((omega*Mu0) / (PI*beta*sum)); }
				tau = Mtau[NLP] = (1.0 / (C*1.0e-12))*(k / beta) * dbdk_bunshi(R, qg, D, w, m, N) / dbdk_bunbo(R, D, w, m, N); // = (1/c)*(dβ/dk) [ps/m]
				modem[NLP] = m;			model[NLP] = l;
				Mbeta[NLP] = beta;

				/* 計算結果の出力 */
				Rinf = R[N] * (bessk(m, w*(Nclad*dr + A)) / bessk(m, w*A));
				fprintf(fp4, "%d, %d, %lf, %lf, %lf, %lf, %lf,", m, l, tau, beta, beta / k, Rinf, eig);
				for (j = 0; j <= N; j++) {
					fprintf(fp4, "%lf,", R2[j]);
					Rlp[NLP][j] = R2[j];
				}

				fprintf(fp4, "\n");
				printf("%d, %d, %lf, %lf, %lf, %lf, %lf\n", m, l, tau, beta, beta / k, eig, R[N]);
				dbeta = k * (n1 - n0) / (double)Nbeta;
				bb = eig;
				l = l + 1;
				NLP = NLP + 1; if (m == 0) { NLP0 = NLP0 + 1; }
				count = 0;
			}

			//  ②「0 < 収束固有値 eig」かつ「-1 < 前回値 bb < 0」であれば
			if (eig > 0.0 && bb < 0.0 && (fabs(bb) < 1.0)) {
				beta = beta + dbeta;
				dbeta = dbeta / 2.0;
				count = count + 1;
			}

			// ③ その他
			else { bb = eig; count = 0; }
			if (count > 3000) {
				dbeta = k * (n1 - n0) / (double)Nbeta;
				beta = beta - dbeta;
			}

		nextin:
			beta = beta - dbeta;
			if (beta < k*n0)	break;
		}
		/****************************************************************************************/
		if (l == 1) { Ntotal = (NLP0) * 2 + (NLP - NLP0) * 4; break; }
	}		// breakとなるのは，mが最高次数に到達したとき
	printf("固有値計算完了\n");

	/* 励振条件設定（MPD配列の算出）*/
	// OFL condition
	if (launch == 0) {
		for (m = 0; m < NLP; m++) {
			if (modem[m] == 0) { MPD[m] = 2.0; }
			else { MPD[m] = 4.0; }
		}
	}	  //縮退数の考慮


		  // RML condition																		
#if 1

	if (launch == 1) {
		double Amev, Amod, Emev, Emod, Ein, xx, rr; //励振条件のとこで出てきた
		int nr;
		for (m = 0; m < NLP; m++) {
			Amev = Amod = 0.0;
			for (i = 0; i < Nxy; i++) {
				xx = (r0 - ((double)Nxy*dx / 2.0)) + ((double)i*dx);
				for (j = 0; j < Nxy; j++) {
					yy = (-((double)Nxy*dx / 2.0)) + ((double)j*dy);
					rr = sqrt(xx*xx + yy * yy);
					nr = (int)(rr / dr);
					if (nr == 0) { Rxy = Rlp[m][0] + (Rlp[m][1] - Rlp[m][0]) * (rr / dr); }	
					else if (nr <= N) { Rxy = Rlp[m][nr] + (Rlp[m][nr + 1] - Rlp[m][nr]) * ((rr - (double)nr*dr) / dr); } //コア
					else if (nr > N) {
						w = aa * sqrt(Mbeta[m] * Mbeta[m] - k * k*n0*n0);
						Rxy = Rlp[m][N] * (bessk(m, w*rr) / bessk(m, w*A));
					}			//クラッド
					Emev = Rxy * cos((double)(modem[m]) * atan2(yy, xx));					
					Emod = Rxy * sin((double)(modem[m]) * atan2(yy, xx));					
					Ein = exp(-((xx - r0)*(xx - r0) + yy * yy) / (w0*w0));					
					Amev = Amev + Emev * Ein*dx*dy;											
					Amod = Amod + Emod * Ein*dx*dy;
				}
			}

			printf("%d ", m);
			Amev = Amev * Amev / (((2.0*omega*Mu0) / Mbeta[m]) * (PI*w0*w0 / 2.0));
			Amod = Amod * Amod / (((2.0*omega*Mu0) / Mbeta[m]) * (PI*w0*w0 / 2.0));			
			MPD[m] = 100.0*(Amev + Amod);
		}
		for (i = 0; i <= NLP - 1; i++) { fprintf(fmpd, "%lf\n", MPD[i]); }
	}
#endif

	printf("%d\n", NLP);


	// RML condition for VCSEL ***********************************************************************
	min = lin = max = nrr1 = nrr2 = 0;
	tauin = betain = Rinfin = eigin = betain_devided_by_k = data = xx1 = xx2 = rr1 = rr2 = 0.0;
	Rinlp = dmatrix(0, (int)(2 * Avin / dr), 0, (int)(2 * Avin / dr));
	MPD2d = dmatrix(0, NLP, 0, NLP);



	modemin = dintvector(0, 10 * Nwkb); init_intvector(modemin, 0, 10 * Nwkb);
	Mbetain = drealvector(0, 10 * Nwkb); init_realvector(Mbetain, 0, 10 * Nwkb);
	
	for (i = 0; i <= (int)(2 * Avin / dr); i++) { for (j = 0; j <= (int)(2 * Avin / dr); j++) { 
		Rinlp[i][j] = 0.0; 		}}
	for (i = 0; i <= NLP; i++) { for (j = 0; j <= NLP; j++) { MPD2d[i][j] = 0.0; } }
	double cef_evev, cef_evod, cef_odev, cef_odod;

	/*********************************************************************************************************************/
	/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
	double Adash;


	if (launch == 2 || launch == 3) {
		//fp9 = fopen("Rinxy.csv", "w");

		if ((fp7 = fopen("FEM_LPmode_VCSELinput.csv", "r")) != NULL) {		//	VCSEL のLPモードの1次元強度分布ファイルを開く

			int minput;
			char MPDfilename[256];
			char coupletype[128];
			if (couple == 0) { sprintf(coupletype, "V2F"); }
			else if (couple == 1) { sprintf(coupletype, "F2F"); }
			else { sprintf(coupletype, "error_warning"); }
			if (couple == 0) { max = num_mode_vin; }
			else if (couple == 1) { max = NLP; }

			for (int ii = 0; ii <= OFFrange / OFFres + 1; ii++) {
				if (launch == 2) { printf("\ncalculating MPD for lateral offset...\n"); }
				if (launch == 3) { printf("\ncalculating MPD for lateral offset of %d um...\n\n", ii * OFFres); }
				fscanf(fp7, "%s", trash);
				double E2m_evin, E2m_ev;
				for (minput = 0; minput < max; minput++) {
					if (ii == 0) {
						fscanf(fp7, "%d, %d, %lf, %lf, %lf, %lf, %lf,", &min, &lin, &tauin, &betain, &betain_devided_by_k, &Rinfin, &eigin);
						modemin[minput] = min;
						modelin[minput] = lin;
						Mbetain[minput] = betain;			}

					printf("%d\t%d\n", modemin[minput], modelin[minput]);
					if (couple == 0) { Adash = Avin; }		if (couple == 1) { Adash = A; }
					
					for (j = 0; j <= (int)(Adash / dr); j++) {
						if (ii == 0) {
							fscanf(fp7, "%lf,", &Rinlp[minput][j]);		}
						//printf("%lf, ", Rinlp[minput][j]);
					}
					if (ii == 0) { fscanf(fp7, "%lf\n", &Rinlp[minput][j]); }
					//printf("%lf\n", Rinlp[minput][j]);


					for (m = 0; m < NLP; m++) {

						Am_evev = 0.0;		Am_evod = 0.0;
						Am_odev = 0.0;		Am_odod = 0.0;
						E2m_evin = 0.0;		E2m_ev = 0.0;
						Rinxy = 0.0;

						for (i = 0; i < Nxy; i++) {
							xx1 = (-(double)Nxy*dx / 2.0) + ((double)i*dx);
							if (launch == 2) {
								xx2 = (r0 - (double)Nxy*dx / 2.0) + ((double)i*dx);
							}
							if (launch == 3) {
								xx2 = ((ii * (double)OFFres * 1.0e-6) - (double)Nxy*dx / 2.0) + ((double)i*dx);
							}

							for (j = 0; j < Nxy; j++) {
								yy = (-(double)Nxy*dx / 2.0) + ((double)j*dy);
								rr1 = sqrt(xx1*xx1 + yy * yy);
								rr2 = sqrt(xx2*xx2 + yy * yy);
								nrr1 = (int)(rr1 / dr);
								nrr2 = (int)(rr2 / dr);
								//printf("%lf\t%d\n", Rinxy, modemin[minput]);

								/*****   ここは入射するレーザ（ファイバ）の，xy平面上の2次元強度分布を計算    *****/
								if (nrr1 == 0) {		// 以下の配列Rinxyは，メモリ内の配列Rlpを使うので読み込みはHDDを介すより速い
									Rinxy = Rinlp[minput][0] + (Rinlp[minput][1] - Rinlp[minput][0]) * (rr1 / dr);
									//printf("%lf\t%d\n", Rinxy, modemin[minput]);
								}
								else if (nrr1 <= Nvin) {
									Rinxy = Rinlp[minput][nrr1] + (Rinlp[minput][nrr1 + 1] - Rinlp[minput][nrr1]) * ((rr1 - (double)nrr1*dr) / dr);
									//printf("%lf\t%d\n", Rinxy, modemin[minput]);
								}
								else if (nrr1 > Nvin) {
									//printf("%d, %lf, %d, %lf\n", minput, Mbetain[minput], Nvin, Avin);
									win = aa * sqrt(Mbetain[minput] * Mbetain[minput] - k * k*n0*n0);
									Rinxy = Rinlp[minput][Nvin] * (bessk(modemin[minput], win*rr1) / bessk(modemin[minput], win*Avin));
									//printf("%e, %e, %e, %e, %e,\t", win, Rinlp[minput][Nvin], bessk(minput, win*rr1)/bessk(minput, win*Avin), rr1, Avin);
									//printf("%lf, %d\n", Rinxy, nrr1);		//printf("%lf\t%d\n", Rinxy, modemin[minput]);
								}
								/*****   ここは受け手側のファイバにおける，xy平面上の2次元強度分布を計算      *****/					/// 以下の配列Rinxyは，メモリ内の配列Rlpを使うので読み込みはHDDを介すより速い
								if (nrr2 == 0) {
									Rxy = Rlp[m][0] + (Rlp[m][1] - Rlp[m][0]) * (rr2 / dr);
								}
								else if (nrr2 <= N) {
									Rxy = Rlp[m][nrr2] + (Rlp[m][nrr2 + 1] - Rlp[m][nrr2]) * ((rr2 - (double)nrr2*dr) / dr);
								}
								else if (nrr2 > N) {
									w = aa * sqrt(Mbeta[m] * Mbeta[m] - k * k*n0*n0);
									Rxy = Rlp[m][N] * (bessk(modem[m], w*rr2) / bessk(modem[m], w*A));
								}

								/*****   重なり積分に用いる電場分布の算出・重なり積分の実行   *****/
								/*if (minput == 1) {
									if (m == 1) {
										printf("%lf, %lf\t%d\n", Rinxy, Rxy, modemin[minput]); }}*/
										//fprintf(fp9, "%lf,", Rinxy);
										/*if (Rinxy == Rxy) { printf("Y "); }
										else printf("N ");*/
								Em_evin = Rinxy * cos((double)(modemin[minput]) * atan2(yy, xx1));
								Em_odin = Rinxy * sin((double)(modemin[minput]) * atan2(yy, xx1));
								Em_ev = Rxy * cos((double)(modem[m]) * atan2(yy, xx2));
								Em_od = Rxy * sin((double)(modem[m]) * atan2(yy, xx2));
								if (i < 1 && j < 1) {
									/*if (modemin[minput] == modem[m]) {
										//printf("%lf, %lf, %lf, %lf\t\t", Em_evin, Em_odin, Em_ev, Em_od); 
										printf("%d   ", m);
										printf("%d, %d, %lf \t %d, %d, %lf\n", modemin[minput], modelin[minput], Mbetain[minput], modem[m], model[m], Mbeta[m]);
									}*/
									if (modem[m] == 0) {
										//printf("%lf, %lf, %lf, %lf\t\t", Em_evin, Em_odin, Em_ev, Em_od); 
										printf("%d ", m);
										//printf("%d, %d, %lf \t %d, %d, %lf\n", modemin[minput], modelin[minput], Mbetain[minput], modem[m], model[m], Mbeta[m]);
									}
									else { printf("%d ", m); }
								}

								E2m_evin = E2m_evin + Em_evin * Em_evin*dx*dy;
								E2m_ev = E2m_ev + Em_ev * Em_ev*dx*dy;
								Am_evev = Am_evev + Em_evin * Em_ev*dx*dy;
								Am_evod = Am_evod + Em_evin * Em_od*dx*dy;
								Am_odev = Am_odev + Em_odin * Em_ev*dx*dy;
								Am_odod = Am_odod + Em_odin * Em_od*dx*dy;
							}

						}

						//printf("%lf\t", E2m_evin);

						cef_evev = Am_evev * Am_evev / (E2m_evin * E2m_ev);
						cef_evod = Am_evod * Am_evod / (E2m_evin * E2m_ev);
						cef_odev = Am_odev * Am_odev / (E2m_evin * E2m_ev);
						cef_odod = Am_odod * Am_odod / (E2m_evin * E2m_ev);

						/*Am_evev = Am_evev * Am_evev / (((2.0*omega*Mu0) / Mbetain[minput]) * (2.0*omega*Mu0) / Mbeta[m]);
						Am_evod = Am_evod * Am_evod / (((2.0*omega*Mu0) / Mbetain[minput]) * (2.0*omega*Mu0) / Mbeta[m]);
						Am_odev = Am_odev * Am_odev / (((2.0*omega*Mu0) / Mbetain[minput]) * (2.0*omega*Mu0) / Mbeta[m]);
						Am_odod = Am_odod * Am_odod / (((2.0*omega*Mu0) / Mbetain[minput]) * (2.0*omega*Mu0) / Mbeta[m]);*/
						if (modem[m] == 0 && modemin[minput] == 0) {
							MPD2d[m][minput] = 100.0*(cef_evev + cef_evod + cef_odev + cef_odod);
						}
						if (modem[m] != 0 && modemin[minput] == 0) {
							MPD2d[m][minput] = 100.0*(cef_evev + cef_evod + cef_odev + cef_odod) / 2;
						}
						if (modem[m] == 0 && modemin[minput] != 0) {
							MPD2d[m][minput] = 100.0*(cef_evev + cef_evod + cef_odev + cef_odod) / 2;
						}
						if (modem[m] != 0 && modemin[minput] != 0) {
							MPD2d[m][minput] = 100.0*(cef_evev + cef_evod + cef_odev + cef_odod) / 2;
						}
						/*if (modem[m] == 0 & modemin[minput]==0) {
							MPD2d[m][minput] = 100.0*(Am_evev + Am_evod + Am_odev + Am_odod);	}
						if (modem[m] != 0 & modemin[minput] == 0){
							MPD2d[m][minput] = 100.0*(Am_evev + Am_evod + Am_odev + Am_odod) / 2;	}
						if (modem[m] == 0 & modemin[minput] != 0){
							MPD2d[m][minput] = 100.0*(Am_evev + Am_evod + Am_odev + Am_odod) / 2;	}
						if (modem[m] != 0 & modemin[minput] != 0) {
							MPD2d[m][minput] = 100.0*(Am_evev + Am_evod + Am_odev + Am_odod) / 4;	}*/

							//printf("%lf, %lf, %lf, %lf\n", Am_evev*100, Am_evod*100, Am_odev*100, Am_odod*100);
					}
					printf("\n");
				}


				printf("\n");
				for (minput = 0; minput < max; minput++) { printf("%d\t%lf\n", modemin[minput], Mbetain[minput]); }
				printf("\n");
				for (m = 0; m < NLP; m++) { printf("%d\t%lf\n", modem[m], Mbeta[m]); }

				if (launch == 2) { sprintf(MPDfilename, "MPD_%s\\MPD2d_%s_%dum_Offset_single(not_dependence_on_x).csv", coupletype, coupletype, r0); }
				if (launch == 3) { sprintf(MPDfilename, "MPD_%s\\MPD2d_%s_%dum_Offset.csv", coupletype, coupletype, ii*OFFres); }


				if ((fp8 = fopen(MPDfilename, "w")) != NULL) {
					fprintf(fp8, "Mode number for receiver-side fiber,");
					if (couple == 0) { fprintf(fp8, "m, l, "); }
					for (j = 0; j < max; j++) { fprintf(fp8, "%d,", j); }
					fprintf(fp8, "\n");
					for (i = 0; i < NLP; i++) {
						fprintf(fp8, "%d,", i);
						if(couple == 0) { fprintf(fp8, "%d, %d, ", modem[i], model[i]); }
						for (j = 0; j < max; j++) {
							// printf("%e, ", MPD2d[i][j]);
							fprintf(fp8, "%e,", MPD2d[i][j]);
						}
						fprintf(fp8, "\n");
					}
					fclose(fp8);
				}
				if (launch == 2) break; 
			}
		}
	}	// if(launch==2) の終了
		/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
		/*********************************************************************************************************************/
	for (m = 0;m < 80;m++) {
		printf("%d, ", modem[m]);
	}



	printf("\nモードパワー分布計算完了\n");

	/** 5. インパルス応答の計算 **/
	taumax = taumin = Mtau[0];
	for (myu = 0; myu < NLP; myu++) {
		if (taumax < Mtau[myu]) { taumax = Mtau[myu]; }
		if (taumin > Mtau[myu]) { taumin = Mtau[myu]; }
	}

	Ntmax = int((taumax*L) / Tv), Ntmin = int((taumin*L) / Tv);
	P = drealvector(0, Ntmax - Ntmin); init_realvector(P, 0, Ntmax - Ntmin); //遅延ステップごとのパワー格納配列
	for (myu = 0; myu < NLP; myu++) {
		x = int((Mtau[myu] * L) / Tv) - Ntmin;
		//P[x] = P[x] +1.0;}
		if (launch == 0) {
			if (modem[myu] == 0) { P[x] = P[x] + 2.0; }
			else { P[x] = P[x] + 4.0; }
		}
		if (launch == 1) {
			//if ( modem[myu] == 1 ) { P[x] = P[x] + 2.0*MPD[myu]; } else { P[x] = P[x] + 4.0*MPD[myu]; }}
			P[x] = P[x] + MPD[myu];
		}
	}

	for (x = 0; x <= Ntmax - Ntmin; x++) { fprintf(fp5, "%lf,%lf\n", double(Ntmin + x)*Tv, P[x]); }
	printf("インパルス応答計算完了\n");

	/** 6. 伝達関数Hwの計算 **/
	/* 周波数応答 M */
	for (y = 0; y <= Nf; y++) {
		Hw = ReHw = ImHw = bw = 0.0;
		for (x = 0; x <= Ntmax - Ntmin; x++) {
			ReHw = ReHw + P[x] * cos((2.0*PI*(fmin + (double)y*dfrq))*(double)x*Tv);
			ImHw = ImHw - P[x] * sin((2.0*PI*(fmin + (double)y*dfrq))*(double)x*Tv);
		}

		M[y] = sqrt(ReHw*ReHw + ImHw * ImHw) / double(Ntotal);
		fprintf(fp6, "%lf,", M[y]);
	}

	/** -3dB帯域 bw **/
	if (M[Nf] < 0.5) {
		for (y = 0; y <= Nf - 1; y++) {
			if (M[y] > 0.5 && M[y + 1] < 0.5) { bw = (fmin*1.0e3) + (dfrq *1.0e3)*(y + (M[y] - 0.5) / (M[y] - M[y + 1])); break; }
		}
	}
	fprintf(fp6, "%lf,", bw);

	printf("周波数応答計算完了\n");

	/** 7. 出射波形の計算 **/
	LT = 2000;
	xmax = Ntmax - Ntmin + 1;
	xL = sizeof(double)*(xmax + LT);
	opulse = drealvector(0, xL); init_realvector(opulse, 0, xL);
	opulse = convolution(xmax, P, LT, pulse);
	for (i = 0; i <= xmax + LT + 1; i++) { fprintf(fopulse, "%lf\n", opulse[i]); }

	/* 配列の記憶領域解放 */
	free_drealvector(GI, 0);
	free_drealvector(q, 0);
	free_drealvector(qg, 0);
	free_drealvector(R, 0);
	free_drealvector(R2, 0);
	free_drealvector(Rb, 0);
	free_drealvector(a, 0);
	free_drealvector(b, 0);
	free_drealvector(ML, 0);
	free_drealvector(MD, 0);
	free_drealvector(Mtau, 0);
	free_dintvector(modem, 0);
	free_drealvector(P, 0);
	free_drealvector(M, 0);
	//
	free_drealvector(Mbeta, 0);
	free_drealvector(MPD, 0);
	free_drealvector(pulse, 0);
	free_drealvector(opulse, 0);

	fprintf(fp2, "Total number of LP modes (WKB), Nwkb,%d\n", Nwkb);
	fprintf(fp2, "Total number of LP modes, NLP,%d\n", NLP);
	fprintf(fp2, "Total number of LP modes with 0th m-order, NLP0,%d\n", NLP0);
	fprintf(fp2, "Total number of modes, Ntotal,%d\n", Ntotal);

	fclose(fp2);
	fclose(fp4);
	fclose(fp5);
	fclose(fp6);	//fclose(fp7);
	fclose(fopulse);
	fclose(fmpd);
	//free_dmatrix(Rlp, 0, sizeofRLP, 0, sizeofRLP);	//free_dmatrix(Rinlp, 0, (int)(2 * A / dr), 0, (int)(2 * A / dr));
	//free_dmatrix(MPD2d, 0, NLP, 0, NLP);

	system("pause");

	return 0;
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
	else printf(" Material number is not correct !\n"); exit(EXIT_FAILURE);					}



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
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch) {
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
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
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
	/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
	free((FREE_ARG)(t[nrl] + ncl - NR_END));
	free((FREE_ARG)(t + nrl - NR_END));
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