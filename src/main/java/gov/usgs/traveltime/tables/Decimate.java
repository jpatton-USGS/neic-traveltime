package gov.usgs.traveltime.tables;

import gov.usgs.traveltime.TauUtil;
import java.util.Arrays;

/**
 * The Decimate class computes a decimation that makes the samples of an array approximately evenly
 * spaced (at a predefined spacing).
 *
 * @author Ray Buland
 */
public class Decimate {
  /** An int containing the current number of data to keep */
  private int numDataToKeep;

  /** An int containing the trial number of data to keep */
  private int trialNumDataToKeep;

  /** A double containing the desired spacing of decimation array values */
  private double targetSpacing;

  /** A double containing the current variance */
  private double currentVariance;

  /** A double[] containing the array to be decimated */
  private double[] decimationArray;

  /**
   * Calculate a decimation for the decimation array such that the differences between the remaining
   * terms is as close to targetSpacing as possible. Note that the first and last elements of
   * decimationArray are always kept. This method figures out the decimation, but doesn't actually
   * implement it. This "slow" method iterates for as long as it takes to minimize the variance
   * between the final grid and the ideal grid.
   *
   * @param decimationArray Array to be decimated
   * @param targetSpacing Target difference between successive elements of decimationArray after
   *     decimation
   * @return An array of booleans, one for each element of decimationArray--if an element is true,
   *     keep the corresponding element of decimationArray
   */
  public boolean[] slowDecimation(double[] decimationArray, double targetSpacing) {
    this.decimationArray = decimationArray;
    this.targetSpacing = targetSpacing;

    // Keep = true if the corresponding decimationArray element will be kept
    boolean[] keep = new boolean[decimationArray.length];
    Arrays.fill(keep, true);

    if (decimationArray.length > 2) {
      // First pass.
      int k0 = 0;
      int k1 = 0;
      int kb = 0;
      double dx1 = 0.0;
      double dx2 = 0.0;
      currentVariance = 0d;
      numDataToKeep = 0;

      for (int j = 1; j < decimationArray.length - 1; j++) {
        dx1 = Math.abs(decimationArray[k1] - decimationArray[j]) - targetSpacing;
        dx2 = Math.abs(decimationArray[k1] - decimationArray[j + 1]) - targetSpacing;

        if (Math.abs(dx2) < Math.abs(dx1)) {
          keep[j] = false;
        } else {
          if (k1 == 0) {
            kb = j;
          }

          k1 = j;
          currentVariance += Math.pow(dx1, 2d);
          numDataToKeep++;
        }
      }
      // Add the last point.
      dx1 =
          Math.abs(decimationArray[k1] - decimationArray[decimationArray.length - 1])
              - targetSpacing;
      currentVariance += Math.pow(dx1, 2d);
      numDataToKeep++;

      if (TablesUtil.deBugLevel > 2) {
        System.out.format("\nInit: %9.3e %9.3e\n", currentVariance / numDataToKeep, doVar(keep));
      }

      // second pass.
      if (numDataToKeep > 1) {
        int pass = 1;
        int nch = 0;

        do {
          k1 = 0;
          int k2 = kb;

          for (int j = kb + 1; j < decimationArray.length; j++) {
            if (keep[j]) {
              k0 = k1;
              k1 = k2;
              k2 = j;
              double var1 = variance(k0, k1, k2, k1 - 1);
              int m1 = trialNumDataToKeep;
              double var2 = variance(k0, k1, k2, k1 + 1);
              int m2 = trialNumDataToKeep;

              if (Math.min(var1 / m1, var2 / m2) < currentVariance / numDataToKeep) {
                // We've reduced the variance.  Decide what to do.
                nch++;
                keep[k1] = !keep[k1];

                // Keep the smallest variance.
                if (var1 / m1 < var2 / m2) {
                  keep[--k1] = true;
                  currentVariance = var1;
                  numDataToKeep = m1;

                  if (TablesUtil.deBugLevel > 2) {
                    System.out.format(
                        "Var1: %9.3e %9.3e %d\n",
                        currentVariance / numDataToKeep, doVar(keep), pass);
                  }
                } else if (var1 / m1 > var2 / m2) {
                  keep[++k1] = true;
                  currentVariance = var2;
                  numDataToKeep = m2;

                  if (TablesUtil.deBugLevel > 2) {
                    System.out.format(
                        "Var2: %9.3e %9.3e %d\n",
                        currentVariance / numDataToKeep, doVar(keep), pass);
                  }
                } else {
                  // If the variances are equal, keep the smallest
                  // number of data.
                  if (m1 <= m2) {
                    keep[--k1] = true;
                    currentVariance = var1;
                    numDataToKeep = m1;

                    if (TablesUtil.deBugLevel > 2) {
                      System.out.format(
                          "M1:   %9.3e %9.3e %d\n",
                          currentVariance / numDataToKeep, doVar(keep), pass);
                    }
                  } else {
                    keep[++k1] = true;
                    currentVariance = var2;
                    numDataToKeep = m2;

                    if (TablesUtil.deBugLevel > 2) {
                      System.out.format(
                          "M2:   %9.3e %9.3e %d\n",
                          currentVariance / numDataToKeep, doVar(keep), pass);
                    }
                  }
                }
              }
              if (k0 == 0) {
                kb = k1;
              }
            }
          }
          pass++;
        } while (nch > 0 && numDataToKeep > 1);
      }
    }

    return keep;
  }

  /**
   * The "fast" method was used in the FORTRAN real-time travel-time calculation for the up-going
   * branches only to save time. It has the advantage of being one pass. There are several
   * differences related to the architecture of the travel-time computation. First, although the
   * decimation is designed to make the spacing in ray travel distance, the distance isn't actually
   * kept and must be estimated from tau. Second, for performance reasons, the algorithm seeks to
   * enforce a minimum distance spacing rather than a uniform target spacing.
   *
   * @param p Normalized ray parameter grid
   * @param tau Normalized tau on the same grid
   * @param xRange Normalized distance at the branch end points
   * @param xMin Normalized minimum distance interval desired
   * @return Decimated, normalized ray parameter grid
   */
  public boolean[] fastDecimation(double[] p, double[] tau, double[] xRange, double xMin) {
    boolean[] keep;
    int n, iBeg, iEnd;
    double xCur, xLast, dx, dx2, sgn, rnd, targetSpacing, xLeast;

    // Scan the current sampling to see if it is already OK.
    xCur = xRange[1];
    for (int i = p.length - 2; i >= 0; i--) {
      xLast = xCur;
      xCur = calcX(p, tau, xRange, i);
      if (Math.abs(xCur - xLast) <= xMin) {

        // It's not OK.  Set up the flag array.
        keep = new boolean[p.length];
        Arrays.fill(keep, true);
        // Set up the decimation algorithm.
        if (Math.abs(xCur - xLast) <= 0.75d * xMin) {
          xCur = xLast;
          i++;
        }
        n = Math.max((int) (Math.abs(xCur - xRange[0]) / xMin + 0.8d), 1);
        dx = (xCur - xRange[0]) / n;
        dx2 = Math.abs(dx / 2d);
        if (dx >= 0d) {
          sgn = 1d;
          rnd = 1d;
        } else {
          sgn = -1d;
          rnd = 0d;
        }
        targetSpacing = xRange[0] + dx;
        iBeg = 1;
        iEnd = 0;
        xLeast = TauUtil.DMAX;

        // Scan the ray parameter grid looking for points to kill.
        for (int j = 1; j <= i; j++) {
          xCur = calcX(p, tau, xRange, j);
          if (sgn * (xCur - targetSpacing) > dx2) {
            // This point looks OK.  See if we have points to kill.
            if (iEnd >= iBeg) {
              for (int k = iBeg; k <= iEnd; k++) keep[k] = false;
            }
            // Reset the kill pointers.
            iBeg = iEnd + 2;
            iEnd = j - 1;
            xLeast = TauUtil.DMAX;
            targetSpacing += (int) ((xCur - targetSpacing - dx2) / dx + rnd) * dx;
          }
          // Look for the best points to kill.
          if (Math.abs(xCur - targetSpacing) < xLeast) {
            xLeast = Math.abs(xCur - targetSpacing);
            iEnd = j - 1;
          }
        }
        // See if there's one more range to kill.
        if (iEnd >= iBeg) {
          for (int k = iBeg; k <= iEnd; k++) keep[k] = false;
        }
        return keep;
      }
    }
    return null;
  }

  /**
   * Return distance as a function of tau (distance is minus the derivative of tau). The method uses
   * a simple three point approximation of the derivative except at the end points where distance is
   * already known.
   *
   * @param p Normalized ray parameter grid
   * @param tau Normalized tau on the same grid
   * @param xRange Normalized distance at the branch end points
   * @param i Grid point where the derivative is required
   * @return Distance corresponding to tau(p_i)
   */
  private double calcX(double[] p, double[] tau, double[] xRange, int i) {
    double h1, h2, hh;

    if (i == 0) return xRange[0];
    else if (i == p.length - 1) return xRange[1];
    else {
      h1 = p[i - 1] - p[i];
      h2 = p[i + 1] - p[i];
      hh = h1 * h2 * (p[i - 1] - p[i + 1]);
      h1 = Math.pow(h1, 2d);
      h2 = -Math.pow(h2, 2d);
      return -(h2 * tau[i - 1] - (h2 + h1) * tau[i] + h1 * tau[i + 1]) / hh;
    }
  }

  /**
   * Compute the variance of various possible values to keep.
   *
   * @param k0 First trial decimation array index
   * @param k1 Second trial decimationArray array index
   * @param k2 Third trial decimationArray array index
   * @param kt Alternate second trial decimationArray array index
   * @return New trial variance of residuals
   */
  private double variance(int k0, int k1, int k2, int kt) {
    double dx1, dx2, newVar;

    dx1 = Math.abs(decimationArray[k0] - decimationArray[k1]) - targetSpacing;
    dx2 = Math.abs(decimationArray[k1] - decimationArray[k2]) - targetSpacing;
    newVar = currentVariance - (Math.pow(dx1, 2d) + Math.pow(dx2, 2d));
    if (kt > k0 && kt < k2) {
      dx1 = Math.abs(decimationArray[k0] - decimationArray[kt]) - targetSpacing;
      dx2 = Math.abs(decimationArray[kt] - decimationArray[k2]) - targetSpacing;
      newVar += Math.pow(dx1, 2d) + Math.pow(dx2, 2d);
      trialNumDataToKeep = numDataToKeep;
    } else {
      dx1 = Math.abs(decimationArray[k0] - decimationArray[k2]) - targetSpacing;
      newVar += Math.pow(dx1, 2d);
      trialNumDataToKeep = numDataToKeep - 1;
    }
    return newVar;
  }

  /**
   * Compute the variance from scratch for testing purposes.
   *
   * @param keep For each element, if true, keep the corresponding decimationArray value
   * @return Variance of absolute differences between kept decimationArray values minus the target
   *     difference
   */
  private double doVar(boolean[] keep) {
    int i = 0, numDataToKeep = 0;
    double currentVariance = 0d;

    for (int j = 1; j < decimationArray.length; j++) {
      if (keep[j]) {
        currentVariance +=
            Math.pow(Math.abs(decimationArray[j] - decimationArray[i]) - targetSpacing, 2d);
        i = j;
        numDataToKeep++;
      }
    }
    return currentVariance / numDataToKeep;
  }
}
