public class Kalman_filter {
    int MP; /* number of measurement vector dimensions */
    int DP; /* number of state vector dimensions */
    int CP; /* number of control vector dimensions */

    public Matrix state_pre; /* predicted state (x'(k)):
                                        x(k)=A*x(k-1)+B*u(k) */
    public Matrix state_post;/* corrected state (x(k)):
                                        x(k)=x'(k)+K(k)*(z(k)-H*x'(k)) */
    public Matrix transition_matrix,transition_matrix_transpose; /* state transition matrix (A) */
    public Matrix control_matrix; /* control matrix (B)
                                       (it is not used if there is no control)*/
    public Matrix control_input; /* control input (u(k))
                                       (it is not used if there is no control)*/
    public Matrix measurement_matrix; /* measurement matrix (H) */
    public Matrix process_noise_cov;  /* process noise covariance matrix (Q) */
    public Matrix measurement_noise_cov; /* measurement noise covariance matrix (R) */
    public Matrix error_cov_pre;  /* priori error estimate covariance matrix (P'(k)):
                                        P'(k)=A*P(k-1)*At + Q)*/
    Matrix gain;  /* Kalman gain matrix (K(k)):
                                        K(k)=P'(k)*Ht*inv(H*P'(k)*Ht+R)*/
    Matrix error_cov_post; /* posteriori error estimate covariance matrix (P(k)):
                                        P(k)=(I-K(k)*H)*P'(k) */
    Matrix temp1; /* temporary matrices */
    Matrix temp2;
    Matrix temp3;
    Matrix temp4;
    Matrix temp5;
    public Kalman_filter() {
          MP = 1;
          DP = 2;
          CP = 0;
          state_pre = new Matrix(DP, 1);
          state_pre.Zero();
          state_post = new Matrix(DP, 1);
          state_post.Zero();
          transition_matrix = new Matrix(DP, DP);
          transition_matrix_transpose = new Matrix(DP, DP);
          transition_matrix.setIdentityMatrix();
          transition_matrix.set(0,1,1);
          process_noise_cov = new Matrix(DP, DP);
          process_noise_cov.setIdentityMatrix();
          measurement_matrix = new Matrix(MP, DP);
//        measurement_matrix.setIdentityMatrix();
          measurement_noise_cov = new Matrix(MP, MP);
          measurement_noise_cov.setIdentityMatrix();
          error_cov_pre = new Matrix(DP, DP);
          error_cov_post = new Matrix(DP, DP);
          error_cov_post.setIdentityMatrix();
          gain = new Matrix(DP, MP);
          if (CP > 0) {
            control_matrix = new Matrix(DP, CP);
            control_matrix.Zero();
          }
    }

    public void predict(){

      /* predicted state (x'(k)):
       x(k)=A*x(k-1)+B*u(k) */
      if(CP > 0){
          state_pre = Matrix.add(Matrix.multiply(transition_matrix,state_post),Matrix.multiply(control_matrix,control_input));
      }else{
          state_pre = Matrix.multiply(transition_matrix,state_post);
      }

      /* P'(k)=A*P(k-1)*At + Q)*/
        temp1 = Matrix.multiply(transition_matrix,error_cov_post);
        Matrix at = Matrix.transpose(transition_matrix);
        error_cov_pre = Matrix.add(Matrix.multiply(temp1,at),process_noise_cov);
//        Matrix result = new Matrix(state_pre);
//        return result;
    }

    private Matrix correct(Matrix measurement){
         /* Kalman gain matrix (K(k)):
        K(k)=P'(k)*Ht*inv(H*P'(k)*Ht+R)*/
        temp2 = Matrix.multiply(error_cov_pre,Matrix.transpose(measurement_matrix));
        temp3 = Matrix.multiply(measurement_matrix,error_cov_pre);
        temp4 = Matrix.add(Matrix.multiply(temp3,Matrix.transpose(measurement_matrix)),measurement_noise_cov);
        gain = Matrix.multiply(temp2,Matrix.invert(temp4));

        /* corrected state (x(k)):
      x(k)=x'(k)+K(k)*(z(k)-H*x'(k)) */
        state_post = Matrix.add(state_pre,
                Matrix.multiply(gain, Matrix.subtract(measurement,Matrix.multiply(measurement_matrix,state_pre))));

        /*  (P(k)):
          P(k)=(I-K(k)*H)*P'(k) */
        error_cov_post = Matrix.subtract(error_cov_pre,
                Matrix.multiply(Matrix.multiply(gain,measurement_matrix),error_cov_pre));

        return state_post;
    }

    public Matrix autoCorrect(Matrix mea){

        predict();
        return correct(mea);

    }


}
