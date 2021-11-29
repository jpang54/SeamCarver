import edu.princeton.cs.algs4.IndexMinPQ;
import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.StdOut;

public class SeamCarver {
    private Picture pic; // defensive copy of input picture

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        if (picture == null) {
            throw new IllegalArgumentException("picture is null");
        }
        pic = new Picture(picture);
    }

    // current picture, give a copy to the client
    public Picture picture() {
        return new Picture(pic);
    }

    // width of current picture
    public int width() {
        return pic.width();
    }

    // height of current picture
    public int height() {
        return pic.height();
    }

    // helper method to calculate x and y gradients, where horizontal and
    // vertical are 0 or 1
    private double calculateEnergy(int col, int row, int horizontal, int vertical) {
        /* @citation Adapted from: https://edstem.org/us/courses/7744/lessons/
        21597/slides/125307. Accessed 11/11/2021. */

        // for the bottom or right pixel, wrap around to 0 if edge is reached
        int rgb1 = pic.getRGB((col + horizontal) % width(),
                              (row + vertical) % height());
        int r1 = (rgb1 >> 16) & 0xFF; // unpack rgb values into r, g, b
        int g1 = (rgb1 >> 8) & 0xFF;
        int b1 = (rgb1) & 0xFF;

        // for the top or left pixel, wrap around to edge-1 if 0 is reached
        int prevCol = col - horizontal;
        int prevRow = row - vertical;
        if (prevCol < 0) {
            prevCol = width() - 1;
        }
        if (prevRow < 0) {
            prevRow = height() - 1;
        }

        int rgb2 = pic.getRGB(prevCol, prevRow);
        int r2 = (rgb2 >> 16) & 0xFF;
        int g2 = (rgb2 >> 8) & 0xFF;
        int b2 = (rgb2) & 0xFF;

        // the given formula for calculating gradients
        return Math.pow(r1 - r2, 2) + Math.pow(g1 - g2, 2) + Math.pow(b1 - b2, 2);

    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (x < 0 || x >= width() || y < 0 || y >= height()) {
            throw new IllegalArgumentException("pixel out of bounds of image");
        }
        double xGradient = calculateEnergy(x, y, 1, 0);
        double yGradient = calculateEnergy(x, y, 0, 1);

        // the given formula for calculating energy
        return Math.sqrt(xGradient + yGradient);
    }

    // helper method to transpose the our stored picture
    private void transpose() {
        Picture transposed = new Picture(height(), width());
        for (int col = 0; col < width(); col++) {
            for (int row = 0; row < height(); row++) {
                transposed.setRGB(row, col, pic.getRGB(col, row));
            }
        }
        pic = transposed;
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        transpose();
        int[] seam = findVerticalSeam();
        transpose();
        return seam;
    }

    // maps from a 2-dimensional (column, row) site to a one-dimensional index
    private int mapper(int col, int row) {
        int oneD = (row * width()) + col;
        return oneD;
    }

    // unmaps from a one-dimensional index to a 2-dimensional (column, row) site
    private int[] unmapper(int oneD) {
        int col = oneD % width();
        int row = (oneD - col) / width();
        int[] coord = new int[2];
        coord[0] = col;
        coord[1] = row;
        return coord;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        /* @citation Adapted from: https://algs4.cs.princeton.edu/44sp/
        DijkstraSP.java.html. Accessed 11/12/2021. */

        // minimum energy to reach a pixel
        double[][] distTo = new double[width()][height()];
        // store the calculated energies for each pixel
        double[][] energy = new double[width()][height()];
        // store the pixels for the minimum energy path
        int[][] vertexTo = new int[width()][height()];
        // the x-coordinates for the pixels in the seam
        int[] seam = new int[height()];

        // initialize all distances to +infinity and uncalculated energies to -1
        for (int i = 0; i < width(); i++) {
            for (int j = 0; j < height(); j++) {
                distTo[i][j] = Double.POSITIVE_INFINITY;
                energy[i][j] = -1;
            }
        }

        // priority queue of pixels
        IndexMinPQ<Double> pq = new IndexMinPQ<Double>(width() * height());

        // insert the first row of pixels and their energy; update arrays
        for (int i = 0; i < width(); i++) {
            energy[i][0] = energy(i, 0);
            distTo[i][0] = energy[i][0];
            pq.insert(i, energy[i][0]);
        }

        // x-coord of pixel in the bottom row on the minimum energy path
        int champIndex = 0;

        // relax pixels in order of their minimum energy path
        while (!pq.isEmpty()) {
            // select unmarked pixel with the smallest energy value
            int coord = pq.delMin();
            int[] unmapped = unmapper(coord); // break 1-D index into 2-D
            int x = unmapped[0];
            int y = unmapped[1];

            // terminate Dijkstra as soon as a pixel on bottom row is deleted
            if (y == height() - 1) {
                champIndex = x;
                break;
            }

            // relax bottom left, directly below, and bottom right of pixel
            if (x - 1 >= 0) {
                relax(x, y, x - 1, pq, distTo, energy, vertexTo);
            }
            relax(x, y, x, pq, distTo, energy, vertexTo);
            if (x + 1 < width()) {
                relax(x, y, x + 1, pq, distTo, energy, vertexTo);
            }
        }

        // retrace the seam backwards
        for (int y = height() - 1; y >= 0; y--) {
            seam[y] = champIndex;
            champIndex = vertexTo[champIndex][y];
        }

        return seam;
    }

    // relax vertex below (xFrom, yFrom) at column xTo and update pq if changed
    private void relax(int xFrom, int yFrom, int xTo,
                       IndexMinPQ<Double> pq, double[][] distTo,
                       double[][] energy, int[][] vertexTo) {
        /* @citation Adapted from: https://algs4.cs.princeton.edu/44sp/
        DijkstraSP.java.html. Accessed 11/12/2021. */

        int yTo = yFrom + 1; // always one row below

        if (energy[xTo][yTo] == -1) {
            energy[xTo][yTo] = energy(xTo, yTo); // store energy if uncalculated
        }

        // if relaxed path is better than current, update arrays and pq
        if (distTo[xTo][yTo] > distTo[xFrom][yFrom] + energy[xTo][yTo]) {
            distTo[xTo][yTo] = distTo[xFrom][yFrom] + energy[xTo][yTo];
            vertexTo[xTo][yTo] = xFrom;
            int oneD = mapper(xTo, yTo);
            if (pq.contains(oneD)) pq.decreaseKey(oneD, distTo[xTo][yTo]);
            else pq.insert(oneD, distTo[xTo][yTo]);
        }
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null) {
            throw new IllegalArgumentException("seam is null");
        }
        if (seam.length != width()) {
            throw new IllegalArgumentException("seam is the wrong length");
        }

        // seam must contain valid y-coordinates with adjacent pixels
        for (int i = 0; i < seam.length; i++) {
            if (seam[i] < 0 || seam[i] >= height()) {
                throw new IllegalArgumentException("array not a valid seam");
            }
            if (i != seam.length - 1) {
                int diff = Math.abs(seam[i + 1] - seam[i]);
                if (diff > 1) {
                    throw new IllegalArgumentException("array not a valid seam");
                }
            }
        }
        if (height() == 1) {
            throw new IllegalArgumentException("height of picture is 1");
        }

        transpose();
        removeVerticalSeam(seam);
        transpose();
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (seam == null) {
            throw new IllegalArgumentException("seam is null");
        }
        if (seam.length != height()) {
            throw new IllegalArgumentException("seam is the wrong length");
        }

        // seam must contain valid x-coordinates with adjacent pixels
        for (int i = 0; i < seam.length; i++) {
            if (seam[i] < 0 || seam[i] >= width()) {
                throw new IllegalArgumentException("array not a valid seam");
            }
            if (i != seam.length - 1) {
                int diff = Math.abs(seam[i + 1] - seam[i]);
                if (diff > 1) {
                    throw new IllegalArgumentException("array not a valid seam");
                }
            }
        }
        if (width() == 1) {
            throw new IllegalArgumentException("width of picture is 1");
        }

        Picture newPic = new Picture(width() - 1, height());

        for (int col = 0; col < width(); col++) {
            for (int row = 0; row < height(); row++) {
                // current pixel on seam, copy pixel to the right
                if (seam[row] == col && col + 1 < width()) {
                    newPic.setRGB(col, row, pic.getRGB(col + 1, row));
                }

                // current pixel not on seam
                if (seam[row] != col && col < width() - 1) {
                    // copy pixel to the right if we passed the seam
                    if (col > seam[row]) {
                        newPic.setRGB(col, row, pic.getRGB(col + 1, row));
                    }
                    else {
                        newPic.setRGB(col, row, pic.getRGB(col, row));
                    }
                }
            }
        }

        pic = newPic;
    }

    public static void main(String[] args) {
        Picture picture1 = new Picture("6x5.png");
        Picture picture2 = new Picture("6x5.png");
        SeamCarver sc1 = new SeamCarver(picture1);
        SeamCarver sc2 = new SeamCarver(picture2);

        StdOut.println(sc1.width()); // 6
        StdOut.println(sc1.height()); // 5
        StdOut.println(sc1.energy(1, 2)); // 138.69

        int[] verticalSeam = sc1.findVerticalSeam();
        StdOut.println("{ ");
        for (int x : verticalSeam)
            StdOut.print(x + " "); // 3 4 3 2 2
        StdOut.println("}");

        int[] horizontalSeam = sc2.findHorizontalSeam();
        StdOut.println("{ ");
        for (int y : horizontalSeam)
            StdOut.print(y + " "); // 2 2 1 2 1 2
        StdOut.println("}");

        sc1.removeVerticalSeam(verticalSeam);
        picture1.show();
        sc1.picture().show();

        sc2.removeHorizontalSeam(horizontalSeam);
        picture2.show();
        sc2.picture().show();
    }
}
