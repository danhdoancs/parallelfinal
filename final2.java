/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package parallelhw5;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import java.util.Arrays;
import org.jfree.ui.RefineryUtilities;

/**
 *
 * @author danh
 */
public class ParallelHw5 {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        // TODO code application logic here

        //PROBLEM 1
        problem1();
        //PROBLEM 2
        //problem2();
        //Problem 3
        //problem3("input2b.txt", "output3.txt");
    }

    public static void problem1() throws IOException {
        GE("input1a.txt", "output1a.txt");
        JOR("input1b.txt", "output1b.txt", false);
    }

    public static void problem2() throws IOException {
        GE("input2a.txt", "output2a.txt");
        JOR("input2b.txt", "output2b.txt", true);
    }

    public static void problem3(String input, String output) throws FileNotFoundException, IOException {

        //a. Calculate the eigen values and vectors
        //init execution time
        //read file
        FileInputStream fileIS = new FileInputStream(input);
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(fileIS));
        String line;

        //get N
        line = bufferedReader.readLine();
        if (!"N".equals(line)) {
            return;
        }
        final int N = Integer.parseInt(bufferedReader.readLine());
        bufferedReader.readLine();

        //get Relaxation parameter Y
        line = bufferedReader.readLine();
        if (!"Y".equals(line)) {
            return;
        }
        String[] stringYs = bufferedReader.readLine().split(" ");
        final float[] Y = new float[stringYs.length];
        for (int i = 0; i < stringYs.length; i++) {
            Y[i] = Float.parseFloat(stringYs[i]);
        }
        bufferedReader.readLine();

        //get Convergence criteria
        line = bufferedReader.readLine();
        if (!"E".equals(line)) {
            return;
        }
        final float E = Float.parseFloat(bufferedReader.readLine());
        bufferedReader.readLine();

        //get Imax
        line = bufferedReader.readLine();
        if (!"Imax".equals(line)) {
            return;
        }
        final int Imax = Integer.parseInt(bufferedReader.readLine());
        bufferedReader.readLine();

        //matrix A
        line = bufferedReader.readLine();
        if (!"A".equals(line)) {
            return;
        }

        final float[][] A = new float[N][N];

        //read each ith row
        for (int i = 0; i < N; i++) {
            String[] row = bufferedReader.readLine().split(" ");
            for (int j = 0; j < N; j++) {
                A[i][j] = Float.parseFloat(row[j]);
            }
        }
        bufferedReader.readLine();

        //get vector b
        line = bufferedReader.readLine();
        if (!"b".equals(line)) {
            return;
        }

        final float[] b = new float[N];
        String[] bArray = bufferedReader.readLine().split(" ");
        for (int j = 0; j < N; j++) {
            b[j] = Float.parseFloat(bArray[j]);
        }

        int[] Inums = new int[Y.length];
        //write output
        File file = new File(output);
        String content = "";

        //init eigen values and vectors
        float[] eigenValues = new float[N];
        float[][] eigenVectors = new float[N][N];
        double bestEigen = 1;
        double smallestEigen = 1;
        double bestY = 1;
        double smallestY = 1;
        boolean isBestYIn2c = false;

        content += "3.a\n";
        for (int y = 0; y < Y.length; y++) {

            float[] oldX = new float[N];
            float[] newX = new float[N];

            //init E
            float[] error = new float[N];

            //initial x(0) values
            for (int i = 0; i < N; i++) {
                oldX[i] = (float) 1;
                error[i] = 1;
            }

            for (Inums[y] = 0; !flagConvergence(error, N, E) && Inums[y] < Imax; Inums[y]++) {
                for (int i = 0; i < N; i++) {
                    int temp = 0;

                    for (int j = 0; j < N; j++) {
                        if (j != i) {
                            temp += oldX[j] * A[i][j];
                        }
                    }

                    newX[i] = (1 - Y[y]) * oldX[i] + (Y[y] / A[i][i]) * (b[i] - temp);

                    error[i] = abs(newX[i] - oldX[i]);
                }

                //newX becomes oldX
                for (int i = 0; i < N; i++) {
                    oldX[i] = newX[i];
                }

                //test
                System.out.println(Arrays.toString(newX));
            }

            //find the eigen values
            float[][] I = new float[N][N];
            float[][] D = new float[N][N];
            float[][] invertedD;
            float[][] B = new float[N][N];
            double[][] M = new double[N][N];
            float[][] Polynomial = new float[N][N];
            float[] lamdas = new float[2];

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (i == j) {
                        I[i][j] = 1;
                        D[i][j] = A[i][j];
                        B[i][j] = 0;
                    } else {
                        I[i][j] = 0;
                        D[i][j] = 0;
                        B[i][j] = A[i][j];
                    }
                }
            }

            invertedD = invert(D);
            float[][] invertedDB = matrixMultiply(invertedD, B, N);

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    M[i][j] = (1 - Y[y]) * I[i][j] + Y[y] * invertedDB[i][j];
                }
            }

            //get eigen values
//            eigenValues[0] = (float) ((M[0][0] + M[1][1]) + sqrt((M[0][0] + M[1][1]) * (M[0][0] + M[1][1]) - 4 * (M[0][0] * M[1][1] - M[0][1] * M[1][0]))) / 2;
//            eigenValues[1] = (float) ((M[0][0] + M[1][1]) - sqrt((M[0][0] + M[1][1]) * (M[0][0] + M[1][1]) - 4 * (M[0][0] * M[1][1] - M[0][1] * M[1][0]))) / 2;
            //get eigen vectors
            // create a symmetric positive definite matrix
            Matrix mM = Matrix.constructWithCopy(M);

            // compute the spectral decomposition
            EigenvalueDecomposition e = mM.eig();
            Matrix mV = e.getV();
            double[] dD = e.getRealEigenvalues();

            //find the best eigen value
            if (dD[0] < bestEigen) {
                bestEigen = dD[0];
                bestY = Y[y];
            }

            content += "Y: " + Y[y] + "\n"
                    + "Eigenvalues: " + vectorToString(dD, N) + "," + "\n"
                    + "Eigenvectors: " + matrixToString(mV.getArray(), N) + "\n";
            ;

            //3.c
            for (int i = 0; i < dD.length; i++) {
                if (dD[i] < smallestEigen) {
                    smallestEigen = dD[i];
                    smallestY = Y[y];
                }
            }
        }

        //print the best Eigen value
        content += "\n3.b\n"
                + "Best eigenvalue: " + bestEigen + " with Y = " + bestY + "\n";

        int minInum = Inums[0];
        double Y2c = 1;
        for (int i = 1; i < Inums.length; i++) {
            if (Inums[i] < minInum) {
                minInum = Inums[i];
                if (bestY == Y[i]) {
                    isBestYIn2c = true;
                    Y2c = Y[i];
                } else {
                    isBestYIn2c = false;
                }
            }
        }

        if (isBestYIn2c) {
            content += "This Y is correspond to the best Y in 2c\n";
        } else {
            content += "This Y is not correspond to the best Y in 2c\n";
        }

        //print the smallest Eigen value
        content += "\n3.c\n"
                + "Smallest eigenvalue: " + smallestEigen + " with Y = " + smallestY + "\n";
        if (Y2c == smallestY) {
            content += "This Y is correspond to the Y in 2c\n";
        } else {
            content += "This Y is not correspond to the Y in 2c\n";
        }

        //3.d
        int Inums2 = 1000;
        float initialValue;
        for(initialValue = -10; initialValue <= 10; initialValue += 0.1) {
            float Yd = (float) smallestY;

            float[] oldX = new float[N];
            float[] newX = new float[N];

            //init E
            float[] error = new float[N];

            //initial x(0) values
            for (int i = 0; i < N; i++) {
                oldX[i] = initialValue;
                error[i] = 1;
            }

            for (Inums2 = 0; !flagConvergence(error, N, E) && Inums2 < Imax; Inums2++) {
                for (int i = 0; i < N; i++) {
                    int temp = 0;

                    for (int j = 0; j < N; j++) {
                        if (j != i) {
                            temp += oldX[j] * A[i][j];
                        }
                    }

                    newX[i] = (1 - Yd) * oldX[i] + (Yd / A[i][i]) * (b[i] - temp);

                    error[i] = abs(newX[i] - oldX[i]);
                }

                //newX becomes oldX
                for (int i = 0; i < N; i++) {
                    oldX[i] = newX[i];
                }
            }
            
            if (Inums2 < minInum)
                break;
        }
        
        content += "\n3.d\n";
        if(Inums2 < minInum)
            content += "Found: x(0) = " + initialValue + " then Inum = " + Inums2;
        else
            content += "Couldn't find.";

        //3.e
        JOR2("input3e.txt", "output3e.txt", true, initialValue);
        
        //create file if not exists
        {
            try (FileOutputStream writer = new FileOutputStream(file)) {
                //create file if not exists
                if (!file.exists()) {
                    file.createNewFile();
                }

                byte[] outputBytes = content.getBytes();

                writer.write(outputBytes);
                writer.flush();
            }
        }
    }

    public static String vectorToString(double[] D, int N) {
        String result = "";
        for (int i = 0; i < N - 1; i++) {
            result += D[i] + ", ";
        }
        result += D[N - 1];
        return result;
    }

    public static String matrixToString(double[][] A, int N) {
        String result = "";
        for (int i = 0; i < N - 1; i++) {
            result += "[";
            for (int j = 0; j < N; j++) {
                result += A[i][j] + ",";
            }
            result += "],";
        }
        result += "[";
        for (int j = 0; j < N; j++) {
            result += A[N - 1][j] + ",";
        }
        result += "]";
        return result;
    }

    public static void JOR(String input, String output, boolean flagPlotChart) throws FileNotFoundException, IOException {

        //read file
        FileInputStream fileIS = new FileInputStream(input);
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(fileIS));
        String line;

        //get N
        line = bufferedReader.readLine();
        if (!"N".equals(line)) {
            return;
        }
        final int N = Integer.parseInt(bufferedReader.readLine());
        bufferedReader.readLine();

        //get Relaxation parameter Y
        line = bufferedReader.readLine();
        if (!"Y".equals(line)) {
            return;
        }
        String[] stringYs = bufferedReader.readLine().split(" ");
        final float[] Y = new float[stringYs.length];
        for (int i = 0; i < stringYs.length; i++) {
            Y[i] = Float.parseFloat(stringYs[i]);
        }
        bufferedReader.readLine();

        //get Convergence criteria
        line = bufferedReader.readLine();
        if (!"E".equals(line)) {
            return;
        }
        final float E = Float.parseFloat(bufferedReader.readLine());
        bufferedReader.readLine();

        //get Imax
        line = bufferedReader.readLine();
        if (!"Imax".equals(line)) {
            return;
        }
        final int Imax = Integer.parseInt(bufferedReader.readLine());
        bufferedReader.readLine();

        //matrix A
        line = bufferedReader.readLine();
        if (!"A".equals(line)) {
            return;
        }

        final float[][] A = new float[N][N];

        //read each ith row
        for (int i = 0; i < N; i++) {
            String[] row = bufferedReader.readLine().split(" ");
            for (int j = 0; j < N; j++) {
                A[i][j] = Float.parseFloat(row[j]);
            }
        }
        bufferedReader.readLine();

        //get vector b
        line = bufferedReader.readLine();
        if (!"b".equals(line)) {
            return;
        }

        final float[] b = new float[N];
        String[] bArray = bufferedReader.readLine().split(" ");
        for (int j = 0; j < N; j++) {
            b[j] = Float.parseFloat(bArray[j]);
        }

        int[] Inums = new int[Y.length];
        //write output
        File file = new File(output);
        String content = "";

        for (int y = 0; y < Y.length; y++) {
            //init execution time
            long startTime = System.currentTimeMillis();

            double[] oldX = new double[N];
            double[] newX = new double[N];

            //init E
            float[] error = new float[N];

            //initial x(0) values
            for (int i = 0; i < N; i++) {
                oldX[i] = (float) 1;
                error[i] = 1;
            }

            for (Inums[y] = 0; !flagConvergence(error, N, E) && Inums[y] < Imax; Inums[y]++) {
                for (int i = 0; i < N; i++) {
                    int temp = 0;

                    for (int j = 0; j < N; j++) {
                        if (j != i) {
                            temp += oldX[j] * A[i][j];
                        }
                    }

                    newX[i] = (1 - Y[y]) * oldX[i] + (Y[y] / A[i][i]) * (b[i] - temp);

                    error[i] = (float)abs(newX[i] - oldX[i]);
                }

                //newX becomes oldX
                for (int i = 0; i < N; i++) {
                    oldX[i] = newX[i];
                }

                //test
                System.out.println(Arrays.toString(newX));
            }

            //stop time
            long stopTime = System.currentTimeMillis();
            long elapsedTime = stopTime - startTime;

            content += "Y = " + Y[y] + "\n"
                    + "Execution time: " + elapsedTime + " miliseconds\n"
                    + "Inum: " + Inums[y] + "\n"
                    + "Solution x: " + Arrays.toString(newX) + "\n";

            if (isTrueSolution(A, newX, b, N)) {
                content += "This solution is right\n";
            } else {
                content += "This solution is wrong\n";
            }
        }

        //create file if not exists
        try (FileOutputStream writer = new FileOutputStream(file)) {
            //create file if not exists
            if (!file.exists()) {
                file.createNewFile();
            }

            byte[] outputBytes = content.getBytes();

            writer.write(outputBytes);
            writer.flush();
        }

        if (flagPlotChart) {
            //create plot
            final BarChart problem2 = new BarChart("Problem 2", Y, Inums);
            problem2.pack();
            RefineryUtilities.centerFrameOnScreen(problem2);
            problem2.setVisible(true);
        }
    }

    public static void JOR2(String input, String output, boolean flagPlotChart, float x0) throws FileNotFoundException, IOException {

        //read file
        FileInputStream fileIS = new FileInputStream(input);
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(fileIS));
        String line;

        //get N
        line = bufferedReader.readLine();
        if (!"N".equals(line)) {
            return;
        }
        final int N = Integer.parseInt(bufferedReader.readLine());
        bufferedReader.readLine();

        //get Relaxation parameter Y
        line = bufferedReader.readLine();
        if (!"Y".equals(line)) {
            return;
        }
        String[] stringYs = bufferedReader.readLine().split(" ");
        final float[] Y = new float[stringYs.length];
        for (int i = 0; i < stringYs.length; i++) {
            Y[i] = Float.parseFloat(stringYs[i]);
        }
        bufferedReader.readLine();

        //get Convergence criteria
        line = bufferedReader.readLine();
        if (!"E".equals(line)) {
            return;
        }
        final float E = Float.parseFloat(bufferedReader.readLine());
        bufferedReader.readLine();

        //get Imax
        line = bufferedReader.readLine();
        if (!"Imax".equals(line)) {
            return;
        }
        final int Imax = Integer.parseInt(bufferedReader.readLine());
        bufferedReader.readLine();

        //matrix A
        line = bufferedReader.readLine();
        if (!"A".equals(line)) {
            return;
        }

        final float[][] A = new float[N][N];

        //read each ith row
        for (int i = 0; i < N; i++) {
            String[] row = bufferedReader.readLine().split(" ");
            for (int j = 0; j < N; j++) {
                A[i][j] = Float.parseFloat(row[j]);
            }
        }
        bufferedReader.readLine();

        //get vector b
        line = bufferedReader.readLine();
        if (!"b".equals(line)) {
            return;
        }

        final float[] b = new float[N];
        String[] bArray = bufferedReader.readLine().split(" ");
        for (int j = 0; j < N; j++) {
            b[j] = Float.parseFloat(bArray[j]);
        }

        int[] Inums = new int[Y.length];
        //write output
        File file = new File(output);
        String content = "";

        for (int y = 0; y < Y.length; y++) {
            //init execution time
            long startTime = System.currentTimeMillis();

            double[] oldX = new double[N];
            double[] newX = new double[N];

            //init E
            float[] error = new float[N];

            //initial x(0) values
            for (int i = 0; i < N; i++) {
                oldX[i] = x0;
                error[i] = 1;
            }

            for (Inums[y] = 0; !flagConvergence(error, N, E) && Inums[y] < Imax; Inums[y]++) {
                for (int i = 0; i < N; i++) {
                    int temp = 0;

                    for (int j = 0; j < N; j++) {
                        if (j != i) {
                            temp += oldX[j] * A[i][j];
                        }
                    }

                    newX[i] = (1 - Y[y]) * oldX[i] + (Y[y] / A[i][i]) * (b[i] - temp);

                    error[i] = (float)abs(newX[i] - oldX[i]);
                }

                //newX becomes oldX
                for (int i = 0; i < N; i++) {
                    oldX[i] = newX[i];
                }

                //test
                System.out.println(Arrays.toString(newX));
            }

            //stop time
            long stopTime = System.currentTimeMillis();
            long elapsedTime = stopTime - startTime;

            content += "Y = " + Y[y] + "\n"
                    + "Execution time: " + elapsedTime + " miliseconds\n"
                    + "Inum: " + Inums[y] + "\n"
                    + "Solution x: " + Arrays.toString(newX) + "\n";

            if (isTrueSolution(A, newX, b, N)) {
                content += "This solution is right\n";
            } else {
                content += "This solution is wrong\n";
            }
        }

        //create file if not exists
        try (FileOutputStream writer = new FileOutputStream(file)) {
            //create file if not exists
            if (!file.exists()) {
                file.createNewFile();
            }

            byte[] outputBytes = content.getBytes();

            writer.write(outputBytes);
            writer.flush();
        }

        if (flagPlotChart) {
            //create plot
            final BarChart problem3 = new BarChart("Problem 3e", Y, Inums);
            problem3.pack();
            RefineryUtilities.centerFrameOnScreen(problem3);
            problem3.setVisible(true);
        }
    }

    public static void GE(String input, String output) throws FileNotFoundException, IOException {

        //init execution time
        long startTime = System.currentTimeMillis();

        //read file
        FileInputStream fileIS = new FileInputStream(input);
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(fileIS));
        String line;

        //get N
        line = bufferedReader.readLine();
        if (!"N".equals(line)) {
            return;
        }
        final int N = Integer.parseInt(bufferedReader.readLine());
        bufferedReader.readLine();

        //matrix A
        line = bufferedReader.readLine();
        if (!"A".equals(line)) {
            return;
        }

        final float[][] A = new float[N][N];
        //matrix Upper
        float[][] U = new float[N][N];
        //matrix Lower
        float[][] L = new float[N][N];

        //read each ith row
        for (int i = 0; i < N; i++) {
            String[] row = bufferedReader.readLine().split(" ");
            for (int j = 0; j < N; j++) {
                A[i][j] = Float.parseFloat(row[j]);

                //copy values to U
                U[i][j] = A[i][j];

                //init values for L
                if (i == j) {
                    L[i][j] = 1;
                } else {
                    L[i][j] = 0;
                }
            }
        }
        bufferedReader.readLine();

        //get vector b
        line = bufferedReader.readLine();
        if (!"b".equals(line)) {
            return;
        }

        final float[] b = new float[N];
        String[] bArray = bufferedReader.readLine().split(" ");
        for (int j = 0; j < N; j++) {
            b[j] = Float.parseFloat(bArray[j]);
        }

        //forward elimination
        for (int k = 0; k < N - 1; k++) {
            for (int i = k + 1; i < N; i++) {
                float multiplier = U[i][k] / U[k][k];

                for (int j = k; j < N; j++) {
                    U[i][j] = U[i][j] - multiplier * U[k][j];
                }

                //find L element values
                L[i][k] = multiplier;

            }

        }

        //solve the equation Ly = b using forward substitution
        float[] y = new float[N];

        for (int i = 0; i < N; i++) {
            y[i] = b[i];

            for (int j = 0; j < i; j++) {
                y[i] -= L[i][j] * y[j];
            }

            y[i] /= L[i][i];
        }

        //solve the equation Ux = y using backward substitution
        double[] x = new double[N];

        for (int i = N - 1; i >= 0; i--) {
            x[i] = y[i];

            for (int j = N - 1; j > i; j--) {
                x[i] -= U[i][j] * x[j];
            }

            x[i] /= U[i][i];
        }

        //stop time
        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;

        //write output
        File file = new File(output);
        String content = "Execution time: " + elapsedTime + " miliseconds\nSolution x: " + Arrays.toString(x);

        if (isTrueSolution(A, x, b, N)) {
            content += "This solution is right";
        } else {
            content += "This solution is wrong";
        }

        //create file if not exists
        try (FileOutputStream writer = new FileOutputStream(file)) {
            //create file if not exists
            if (!file.exists()) {
                file.createNewFile();
            }

            byte[] outputBytes = content.getBytes();

            writer.write(outputBytes);
            writer.flush();
        }

    }

    public static boolean flagConvergence(float[] errors, int errorSize, float E) {
        for (int i = 0; i < errorSize; i++) {
            if (errors[i] > E) {
                return false;
            }
        }
        return true;
    }

    public static boolean isTrueSolution(float[][] A, double[] x, float[] B, int N) {

        for (int i = 0; i < N; i++) {
            float answer = 0;
            for (int j = 0; j < N; j++) {
                answer += A[i][j] * x[j];
            }

            if (abs(answer - B[i]) > 1) {
                return false;
            }
        }

        return true;
    }

    public static float[][] invert(float a[][]) {

        int n = a.length;

        float x[][] = new float[n][n];

        float b[][] = new float[n][n];

        int index[] = new int[n];

        for (int i = 0; i < n; ++i) {
            b[i][i] = 1;
        }

        // Transform the matrix into an upper triangle
        gaussian(a, index);

        // Update the matrix b[i][j] with the ratios stored
        for (int i = 0; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    b[index[j]][k]
                            -= a[index[j]][i] * b[index[i]][k];
                }
            }
        }

        // Perform backward substitutions
        for (int i = 0; i < n; ++i) {

            x[n - 1][i] = b[index[n - 1]][i] / a[index[n - 1]][n - 1];

            for (int j = n - 2; j >= 0; --j) {

                x[j][i] = b[index[j]][i];

                for (int k = j + 1; k < n; ++k) {

                    x[j][i] -= a[index[j]][k] * x[k][i];

                }

                x[j][i] /= a[index[j]][j];

            }

        }

        return x;

    }

    public static void gaussian(float a[][], int index[]) {

        int n = index.length;

        float c[] = new float[n];

        // Initialize the index
        for (int i = 0; i < n; ++i) {
            index[i] = i;
        }

        // Find the rescaling factors, one from each row
        for (int i = 0; i < n; ++i) {

            float c1 = 0;

            for (int j = 0; j < n; ++j) {

                float c0 = Math.abs(a[i][j]);

                if (c0 > c1) {
                    c1 = c0;
                }

            }

            c[i] = c1;

        }

        // Search the pivoting element from each column
        int k = 0;

        for (int j = 0; j < n - 1; ++j) {

            float pi1 = 0;

            for (int i = j; i < n; ++i) {

                float pi0 = Math.abs(a[index[i]][j]);

                pi0 /= c[index[i]];

                if (pi0 > pi1) {

                    pi1 = pi0;

                    k = i;

                }

            }

            // Interchange rows according to the pivoting order
            int itmp = index[j];

            index[j] = index[k];

            index[k] = itmp;

            for (int i = j + 1; i < n; ++i) {

                float pj = a[index[i]][j] / a[index[j]][j];

                // Record pivoting ratios below the diagonal
                a[index[i]][j] = pj;

                // Modify other elements accordingly
                for (int l = j + 1; l < n; ++l) {
                    a[index[i]][l] -= pj * a[index[j]][l];
                }

            }

        }

    }

    public static float determinant(float A[][], int N) {

        float det = 0;

        if (N == 1) {

            det = A[0][0];

        } else if (N == 2) {

            det = A[0][0] * A[1][1] - A[1][0] * A[0][1];

        } else {

            det = 0;

            for (int j1 = 0; j1 < N; j1++) {

                float[][] m = new float[N - 1][];

                for (int k = 0; k < (N - 1); k++) {

                    m[k] = new float[N - 1];

                }

                for (int i = 1; i < N; i++) {

                    int j2 = 0;

                    for (int j = 0; j < N; j++) {

                        if (j == j1) {
                            continue;
                        }

                        m[i - 1][j2] = A[i][j];

                        j2++;

                    }

                }

                det += Math.pow(-1.0, 1.0 + j1 + 1.0) * A[0][j1] * determinant(m, N - 1);

            }

        }

        return det;

    }

    public static float[][] matrixMultiply(float[][] A, float[][] B, int N) {
        float[][] AB = new float[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    AB[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return AB;
    }
}
