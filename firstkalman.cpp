#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

//[[Rcpp::export]]
class Kalman {
	private:
		int t;     //time index
		mat X;  //state vector
		mat P;  //covariance matrix of state vector
		mat predX;
		mat predP;
		mat PHI; //state transition matrix
		mat UPSILON;
		mat GAMMA;
		mat Q; //covariance of state vector
		mat R; //covariance of observation vector
		mat GAIN;
	public:	
		Kalman(mat mu0, mat Sigma0, mat Phi0, mat Q0, mat R0); //constructor
		void updateState(mat u, mat y, mat A);  //update state vector once every time period
		mat pointPredict(mat aheadA);
		mat predError(mat aheadA);
		//maybe include some methods to output stuff
};

Kalman::Kalman(mat mu0, mat Sigma0, mat Phi0, mat Q0, mat R0){
	X = mu0;
	P = Sigma0;
	Q = Q0;
	R = R0;
	PHI = Phi0;
	t = 0;
}

//make sure all the As are the same
void Kalman::updateState(mat u, mat y, mat A){  //doesn't store y or u
	int p = P.n_rows;
  mat id = eye(p,p);
	predX = PHI*X + UPSILON*u;               //make sure you only call once
	predP = PHI*P*(trans(PHI)) + Q;
	mat GAIN = predP*trans(A)*(A*predP*trans(A)+R).i();
	X = predX + GAIN*(y - A*predX - GAMMA*u);
	P = (id - GAIN*A)*predP;
	t += 1;
}

//make sure you don't screw up the indices
mat Kalman::pointPredict(mat aheadA){
	mat predY = aheadA*predX;
	return predY;
}

mat Kalman::predError(mat aheadA){
	mat predVarMat = aheadA*predP*trans(aheadA) + R;
	return predVarMat;
}

//there are probably going to be some issues abotu the order you can call things in
