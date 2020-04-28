import org.omg.CORBA.MARSHAL;

import java.util.*;

public class test {
    public static void main(String[] args) {

        Kalman_filter kalman_filter = new Kalman_filter();
        Random rand = new Random();

        kalman_filter.state_post.set(0,0,100);
        kalman_filter.state_post.set(1,0,10);
        kalman_filter.measurement_matrix.setAll(1,0);
        kalman_filter.process_noise_cov.scale(1e-5);
        kalman_filter.measurement_noise_cov.scale(1e-1);

        for(int x = 0; x<100;x++){
            kalman_filter.transition_matrix.setAll(1,x,0,1);
            double y = 100 +10*x;
            double noise = 10 * Math.sin((40.0 * Math.PI / (float)200) * x);
            Matrix z_k = new Matrix(1,1);
            z_k.setAll(y + noise);
            System.out.println(kalman_filter.autoCorrect(z_k).get(0,0)+" "+z_k.get(0,0));
        }



    }

}
