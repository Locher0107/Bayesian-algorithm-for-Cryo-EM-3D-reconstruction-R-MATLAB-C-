


#include <iostream>  
#include <armadillo>  /*you have to download it*/
#include<chrono>
#include<random>
#include<boost\math\distributions\gamma.hpp>
#include<fstream>
#include<string>

using namespace std;
using namespace arma;
using namespace boost::math;


//generate student-t random variable
//with df=df_0 and length ns
//t-distribution is not related to our algorithm
vec rstudent(int ns, double df_0, double scale) {
	auto seed_t = chrono::high_resolution_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed_t);
	std::student_t_distribution<double> t_distribution(df_0);

	vec t_store(ns);
	for (int i = 0; i < ns; i++) {
		t_store(i) = t_distribution(generator)*scale;
	}
	return t_store;
}

//generate gamma random variable
double rgamma(double shape, double scale) {

	auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
	std::tr1::mt19937 eng(seed);
	std::tr1::gamma_distribution<double> gamma(shape);
	return scale*gamma(eng);
}

//draw beta and generate multivariate Normal distribution
vec draw_beta(colvec mu_n, mat Lambda_n_inv, double sigma) {
	//get the length of vector;
	arma::arma_rng::set_seed_random();
	int ncol = mu_n.n_elem;
	//Lambda_n_inv.print("Lambda_n_inv:");
	mat L0 = chol(pow(sigma, 2)*Lambda_n_inv);
	//L0.print("L0:");
	//generate N(0,1)
	colvec X0 = randn<vec>(ncol);
	//X0.print("X0:");
	vec beta0 = L0*X0 + mu_n;
	//beta0.print("beta0:");
	return beta0;
}


vec draw_lambda(vec lambda, vec beta, double sigma, double tau) {
	arma::arma_rng::set_seed_random();
	vec gamma_0 = 1 / pow(lambda, 2);
	vec mu_hat = beta / (tau*sigma);
	vec mu = randu<vec>(lambda.n_elem);
	mu = mu % (1 / (1 + gamma_0));
	//mu.print("mu");
	/*generate gamma with non-standard distribution*/
	/*get cdf first*/
	
	vec upper_gamma = 1 - exp(  -1.0 / 2 * square(mu_hat)% ((1 - mu) / mu)); //???
	//mu_hat.print("mu_hat");
	/*F(X) is uniformly distributed*/
	//upper_gamma.print("upper_gamma");
	vec mu2 = randu(lambda.n_elem) % upper_gamma;
	//mu2.print("mu2");
	/*find the inverse quantile*/
	vec gamma_final = log(1 - mu2) / (-1.0 / 2 * square(mu_hat));
	//gamma_final.print("gamma_final");
	return 1 / sqrt(gamma_final);

}


double draw_tau(vec lambda, vec beta, double sigma, double tau) {
	arma::arma_rng::set_seed_random();
	double gamma_0 = 1.0 / pow(tau,2);
	double mu_hat_2 = sum(square((beta*sigma)/lambda));
	vec t0 = randu(1);
	double mu = (t0*1.0 / (1 + gamma_0))[0];
	int p = lambda.n_elem;
	double rate_tau = mu_hat_2 / 2.0;
	/*generate gamma with non-standard distribution*/
	double upper_gamma = gamma_p((p + 1.0) / 2, rate_tau*(1 - mu) / mu);
	double mu_2 = as_scalar(randu(1))*upper_gamma;
	double gamma_final = gamma_p_inv((p + 1.0) / 2, mu_2)*(1 / rate_tau);
	return 1 / sqrt(gamma_final);

}

double draw_sigma(double a_n, double b_0, double yty, colvec mu_n, mat Lambda_n) {
	//mat z_n =trans(mu_n)*Lambda_n*mu_n;
	//double z_n0 = z_n[0][0];
	double b_n = b_0 + (0.5*(yty - as_scalar(trans(mu_n)*Lambda_n*mu_n)));
	/*vec value0= randg<vec>(1, distr_param(a_n, 1 / b_n)) ;*/

	return sqrt(1 / rgamma(a_n, 1.0 / b_n));

}

//make horseshoe as a tuple. So we can
//return(mat beta, mat lambda, vec tau, vec sigma)
std::tuple<mat, mat, vec, vec> horseshoe_main(vec y, mat x, int iter = 2000, vec ab = { 1.0, 1.0 }) {

	arma::arma_rng::set_seed_random();
	/*parameter storage*/
	mat beta_bayes_hs(iter , x.n_cols);
	mat lambda(iter+10 , x.n_cols);
	vec sigma(iter+10);
	vec tau(iter+10 );
	mat Lambda0(x.n_cols, x.n_cols);


	double a_0 = ab(0);
	double b_0 = ab(1);
	/*initial value*/
	beta_bayes_hs.row(0).fill(0);
	lambda.row(0).fill(1);
	tau(0) = 1.0;
	sigma(0) = 1.0;
	Lambda0.fill(0);

	//precaculatiyton
	mat XtX = x.t()*x;
	vec Xty = x.t()*y;
	double yty = as_scalar(y.t()*y);
	double a_n = a_0 + y.n_elem / 2;


	//sampling and iteration
	for (int it = 1; it < iter; it++) {
		//know the index of iteration
		//cout << "iteration:   " << it << endl;

		//To get Lambda0, Lambda_n, Lambda_n_inv,mu_n
		Lambda0.diag() = 1 / pow(lambda.row(it - 1)*tau(it - 1), 2);
		mat Lambda_n = XtX + Lambda0;
		mat Lambda_n_inv = inv(Lambda_n);
		vec mu_n = Lambda_n_inv*Xty;


		/*******************************
		*********sampling***************
		********************************/
		//sampling beta
		vec draw_beta_store = draw_beta(mu_n, Lambda_n_inv, sigma(it - 1));
		//because draw_beta_store is column vector, 
		//beta_bayes_hs.row(it) is row vector, we need to 
		//take the transfer of column vector
		beta_bayes_hs.row(it) = draw_beta_store.t();

		//sampling sigma
		sigma(it) = draw_sigma(a_n, b_0, yty, mu_n, Lambda_n);

		//sampling lambda
		 vec draw_lambda_store= draw_lambda(lambda.row(it - 1).t(), beta_bayes_hs.row(it).t(), sigma(it), tau(it - 1));
		// draw_lambda_store.print();
		 lambda.row(it) = draw_lambda_store.t();

		//sampling tau
		tau(it) = draw_tau(lambda.row(it).t(), beta_bayes_hs.row(it).t(), sigma(it), tau(it - 1));

	}

	//return 
	return std::make_tuple(beta_bayes_hs, lambda, tau, sigma);



}

/*std::tuple<int, int> foo_tuple()
{
	// Error until C++17
	return std::make_tuple(1, -1); // Always works
}*/

int main()
{	
	// Record start time
	auto start = std::chrono::high_resolution_clock::now();

	//vec t0 = randu(1);
	//vec b = t0 * 1;
	//b.print();

   //test draw function
/*	//dar_beta checked!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	vec mu_n = { 1,1 };
	mat Lambda_n_inv(2, 2); Lambda_n_inv.eye();
	double sigma = 1;
	vec btee_ini = { 0,0 };
	mat whole(2, 10000);
	for (int i = 0; i < 10000; i++) {
		whole.col(i) = draw_beta(mu_n, Lambda_n_inv, sigma);
		//btee_ini = btee_ini + btee;
	}
	cout << "Var" << var(whole.row(1)) << endl;
*/
/*	//test rgamma
	vec whole(1000000);
	for (int i = 0; i < 1000000; i++) {
		whole(i) = rgamma(1.0,1.0);
		//btee_ini = btee_ini + btee;
	}
	cout << "Mean" << mean(whole) << endl;
*/
	ofstream outFile;
	outFile.open("sample_05_200size.txt");


	vec mle(500); vec horseshoe_sse(500);
	for (int dat_it = 0; dat_it < 1; dat_it++) {
		//span matrix
		/*
		mat attt=randu<mat>(5,5);
		attt.print("attt");
		attt(span(2, 4), span(2, 4)).print("span matrix");
		*/

		//main algorithm
		//generate simulation data
		arma_rng::set_seed_random();

		cout << "iteration:  " << dat_it+1 << endl;
		outFile << "iteration:  " << dat_it+1 << endl;


		double w = 0.2;
		int ns = 200;
		int ksi = 2;
		vec ones_(200); ones_.fill(1.0);
		vec prop = floor(ones_-randu<vec>(ns) +w);
		vec theta = prop%rstudent(ns, ksi, 3);
		/*for (int i = 0; i < 200; i++) {
			cout <<"theta   "<< theta(i) << endl;
		}*/
	
	
		//theta.print("Theta =\n");
		//generate X_data
		mat X_data = diagmat(ones<vec>(ns));
		vec y = vectorise(X_data*theta) + randn<vec>(ns);

		//******************************************
		//***********************1. test MLE
		//*****************************************
		vec beta_hat = inv(X_data.t()*X_data)*X_data.t()*y;
		mle(dat_it) = sum(square(beta_hat - theta));
		cout << "The sum of square error of MLE is :\n" << mle(dat_it) << endl;
		outFile << "The sum of square error of MLE is :\n" << mle(dat_it) << endl;

		//*****************************88
		//*****************TEST Horseshoe
		//***********************************
		mat beta_bayes_hs_T; mat lambda_T;
		vec tau_T; vec sigma_T;
		std::tie(beta_bayes_hs_T, lambda_T, tau_T, sigma_T) = horseshoe_main(y, X_data, 8000);
		//dim=0 means operates on all elements in a column
		//burn-in start at 1001
		mat beta_hs = mean(beta_bayes_hs_T(span(1000, 1999), span::all), 0);
		//cout << arma::size(beta_hs) << endl;

		horseshoe_sse(dat_it) = sum(square(beta_hs.t() - theta));
		cout << "The sum of square error of Horseshoe is :\n" << horseshoe_sse(dat_it) << endl;
		outFile << "The sum of square error of Horseshoe is :\n" << horseshoe_sse(dat_it) << endl;
	}

	/*for (int i = 0; i < 5; i++) {
		cout <<"mle"<< mle(i) << endl;
		cout <<"horseshoe"<< horseshoe_sse(i) << endl;
	}*/
	
	double mle_mse = mean(mle(span(0,49)));
	double horseshoe_mse = mean(horseshoe_sse(span(0,49)  ));
	cout << "Average of SSE for MLE: " << mle_mse << endl;
	cout << "Average of SSE for Horseshoe: " << horseshoe_mse << endl;
	outFile << "Average of SSE for MLE: " << mle_mse << endl;
	outFile<< "Average of SSE for Horseshoe: " << horseshoe_mse << endl;

	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s\n";
	outFile << "Elapsed time: " << elapsed.count() << " s\n";

	while (1);
	return 0;
}



	//arma::mat X;
	//arma::vec beta;

	//beta.resize(2);

	//beta(0) = 1.0;
	//beta(1) = 3.0;

	//X.resize(3, 2);

	//X(0, 0) = 1.0;
	//X(0, 1) = 2.0;
	//X(1, 0) = 3.0;
	//X(1, 1) = 4.0;
	//X(2, 0) = 5.0;
	//X(2, 1) = 6.0;

	//X.print();
	//mat colx(2, 1);
	//colx(0, 0) = 1;
	//colx(1, 0) = 1;
	//colx.print();

	//beta.print();
	//double sigma_test = 3;
	//mat beta_mat(beta);
	//mat xcolx = X*colx;
	//beta_mat.print();
	////mat XBETA = X*beta_mat;
	////std::cout << X * beta_mat << std::endl;



	//mat A = randu<mat>(2, 2);
	//mat C = randu<mat>(2, 2);
	//mat U = A * C;
	////mat pld = chol(A);
	//A.print("A:");
	////pld.print("pld:");
	//mat V = inv(U);
	//V.print("V:");
	//cout << rgamma(5, 3) << endl;

	////cholesky decomposition
	//mat BBB = randu<mat>(5, 5);
	//mat yyy = BBB.t()*BBB;
	//yyy.print("yyy:");

	//rowvec vec_test = { 1,2,3 };
	//vec_test.print("vec_test:");


	//mat R1 = chol(yyy);
	//R1.print("R1=");
	////mat R2 = chol(yyy, "lower");

	////***************************
	////*******test Draw Function********//////////
	////**************************************/////
	////test draw_beta
	///*vec mu_N = { 2,1 };
	////generate Matrix
	//mat Sigma_0(2, 2);
	//Sigma_0(0,0)=1;
	//Sigma_0(0, 1) = 0;
	//Sigma_0(1, 0) = 0;
	//Sigma_0(1, 1) = 1;
	//Sigma_0.print("Sigma_0=");
	//vec accum = {0,0};
	////to caculate the mean
	//for (int i = 0; i<1000; i++) {
	//	accum = accum + draw_beta(mu_N, Sigma_0, 1);
	//}
	//cout << "mean\n" << accum / 1000 << endl;

	////*****************************************
	////*************************************




	///***********************************************
	//******************assign value test***********
	//**********************************************/

	//mat betabayes(2, 2);
	//betabayes.row(0).fill(1);
	//betabayes.row(1).fill(1);
	//betabayes.print("betabayes:");
	////betabayes.row(1)= 0;
	//rowvec r(5); r.fill(123.0);
	//r.print("r:");
	////betabayes.print("betabeyes:");


	//vec lambda_vec = { 1,1,2 };
	//double tau_0 = 2;
	//(1 / pow(lambda_vec*tau_0, 2)).print("vec: ");

	////diaganol assign
	//mat diag_A = randu<mat>(5, 5);
	//diag_A.print("diag_A_origin");
	//rowvec diag_a(5); diag_a.fill(1);
	//diag_A.diag() = diag_a;
	//diag_A.print("print again:");


	////test gamma cdf function
	//double scale_0 = double(1.0 / 5);
	//cout << "scale_0:" << scale_0 << endl;
	//cout << "gamma cdf:\n" << gamma_p(1, scale_0) << endl;
	//double rrra = as_scalar(randu(1));
	//cout << "uniform distribution\n" << rrra << endl;





/*
std::tuple<double, char, std::string> get_student(int id)
{
	if (id == 0) return std::make_tuple(3.8, 'A', "Lisa Simpson");
	if (id == 1) return std::make_tuple(2.9, 'C', "Milhouse Van Houten");
	if (id == 2) return std::make_tuple(1.7, 'D', "Ralph Wiggum");
	throw std::invalid_argument("id");
}

int main()
{
	auto student0 = get_student(0);
	std::cout << "ID: 0, "
		<< "GPA: " << std::get<0>(student0) << ", "
		<< "grade: " << std::get<1>(student0) << ", "
		<< "name: " << std::get<2>(student0) << '\n';

	double gpa1;
	char grade1;
	std::string name1;
	std::tie(gpa1, grade1, name1) = get_student(1);
	std::cout << "ID: 1, "
		<< "GPA: " << gpa1 << ", "
		<< "grade: " << grade1 << ", "
		<< "name: " << name1 << '\n';


	while (1);
	return 0;
}*/