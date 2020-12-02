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
  /**
   * A integer value containing the current number of data to keep. This is global to the class
   * because of the variance() function
   */
  private int currentNumToKeep;

  /**
   * A integer value containing the trial number of data to keep. This is global to the class
   * because of the variance() function
   */
  private int trialNumToKeep;

  /**
   * A Double value containing the desired spacing of array values. This is global to the class
   * because of the variance() function
   */
  private double targetSpacing;

  /**
   * A Double value containing the current variance of residuals. This is global to the class
   * because of the variance() function
   */
  private double currentVariance;

  /** A Double array containing the array of values to be decimated */
  private double[] decimationArray;

  /**
   * This method calculates a decimation for the array decimationArray such that the differences
   * between the remaining terms is as close to targetSpacing as possible. Note that the first and
   * last elements of decimationArray are always kept. This method figures out the decimation, but
   * doesn't actually implement it. This "slow" method iterates for as long as it takes to minimize
   * the variance between the final grid and the ideal grid.
   *
   * @param decimationArray Array to be decimated
   * @param targetSpacing A double value containing the target difference between successive
   *     elements of decimationArray after decimation
   * @return An array of booleans, one for each element of decimationArray, if an element is true,
   *     keep the corresponding element of decimationArray
   */
  public boolean[] slowDecimation(double[] decimationArray, double targetSpacing) {
    int k0, k1, k2, kb = 0, nch, m1, m2, pass;
    double dx1, dx2, var1, var2;

    this.decimationArray = decimationArray;
    this.targetSpacing = targetSpacing;

    boolean[] keep = new boolean[decimationArray.length];
    // True if this decimationArray element will be kept
    Arrays.fill(keep, true);

    if (decimationArray.length > 2) {
      // First pass.
      k1 = 0;
      currentVariance = 0d;
      currentNumToKeep = 0;

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
          currentNumToKeep++;
        }
      }

      // Add the last point.
      dx1 =
          Math.abs(decimationArray[k1] - decimationArray[decimationArray.length - 1])
              - targetSpacing;
      currentVariance += Math.pow(dx1, 2d);
      currentNumToKeep++;

      if (TablesUtil.deBugLevel > 2) {
        System.out.format("\nInit: %9.3e %9.3e\n", currentVariance / currentNumToKeep, doVar(keep));
      }

      // second pass.
      if (currentNumToKeep > 1) {
        pass = 1;
        do {
          k1 = 0;
          k2 = kb;
          nch = 0;
          for (int j = kb + 1; j < decimationArray.length; j++) {
            if (keep[j]) {
              k0 = k1;
              k1 = k2;
              k2 = j;
              var1 = variance(k0, k1, k2, k1 - 1);
              m1 = trialNumToKeep;
              var2 = variance(k0, k1, k2, k1 + 1);
              m2 = trialNumToKeep;

              if (Math.min(var1 / m1, var2 / m2) < currentVariance / currentNumToKeep) {
                // We've reduced the variance.  Decide what to do.
                nch++;
                keep[k1] = !keep[k1];

                // Keep the smallest variance.
                if (var1 / m1 < var2 / m2) {
                  keep[--k1] = true;
                  currentVariance = var1;
                  currentNumToKeep = m1;

                  if (TablesUtil.deBugLevel > 2) {
                    System.out.format(
                        "Var1: %9.3e %9.3e %d\n",
                        currentVariance / currentNumToKeep, doVar(keep), pass);
                  }
                } else if (var1 / m1 > var2 / m2) {
                  keep[++k1] = true;
                  currentVariance = var2;
                  currentNumToKeep = m2;

                  if (TablesUtil.deBugLevel > 2) {
                    System.out.format(
                        "Var2: %9.3e %9.3e %d\n",
                        currentVariance / currentNumToKeep, doVar(keep), pass);
                  }
                } else {
                  // If the variances are equal, keep the smallest
                  // number of data.
                  if (m1 <= m2) {
                    keep[--k1] = true;
                    currentVariance = var1;
                    currentNumToKeep = m1;

                    if (TablesUtil.deBugLevel > 2) {
                      System.out.format(
                          "M1:   %9.3e %9.3e %d\n",
                          currentVariance / currentNumToKeep, doVar(keep), pass);
                    }
                  } else {
                    keep[++k1] = true;
                    currentVariance = var2;
                    currentNumToKeep = m2;

                    if (TablesUtil.deBugLevel > 2) {
                      System.out.format(
                          "M2:   %9.3e %9.3e %d\n",
                          currentVariance / currentNumToKeep, doVar(keep), pass);
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
        } while (nch > 0 && currentNumToKeep > 1);
      }
    }

    return keep;
  }

  /**
   * This method calculates a decimation for the array rayParamGrid using the "fast" method
   *
   * <p>The "fast" method was used in the FORTRAN real-time travel-time calculation for the up-going
   * branches only to save time. It has the advantage of being one pass. There are several
   * differences related to the architecture of the travel-time computation. First, although the
   * decimation is designed to make the spacing in ray travel distance, the distance isn't actually
   * kept and must be estimated from tau. Second, for performance reasons, the algorithm seeks to
   * enforce a minimum distance spacing rather than a uniform target spacing.
   *
   * @param rayParamGrid An array of doubles containing the normalized ray parameter grid
   * @param tauGrid An array of doubles containing the normalized tau on the same grid
   * @param distanceRange An array of doubles containing the normalized distance at the branch end
   *     points. Note that the length of this array is expected to be two, with the 0 index
   *     containing the minumum range point, and 1 containing the maximum range point.
   * @param distanceMin A double value contining the normalized minimum distance interval desired
   * @return An array of booleans containing the decimated, normalized ray parameter grid
   */
  public boolean[] fastDecimation(
      double[] rayParamGrid, double[] tauGrid, double[] distanceRange, double distanceMin) {
    boolean[] keep;
    int n, iBeg, iEnd;
    double xCur, xLast, dx, dx2, sgn, rnd, targetSpacing, xLeast;

    // Scan the current sampling to see if it is already OK.
    xCur = distanceRange[1];
    for (int i = rayParamGrid.length - 2; i >= 0; i--) {
      xLast = xCur;
      xCur = calcX(rayParamGrid, tauGrid, distanceRange, i);

      if (Math.abs(xCur - xLast) <= distanceMin) {
        // It's not OK.  Set up the flag array.
        keep = new boolean[rayParamGrid.length];
        Arrays.fill(keep, true);

        // Set up the decimation algorithm.
        if (Math.abs(xCur - xLast) <= 0.75d * distanceMin) {
          xCur = xLast;
          i++;
        }

        n = Math.max((int) (Math.abs(xCur - distanceRange[0]) / distanceMin + 0.8d), 1);
        dx = (xCur - distanceRange[0]) / n;
        dx2 = Math.abs(dx / 2d);

        if (dx >= 0d) {
          sgn = 1d;
          rnd = 1d;
        } else {
          sgn = -1d;
          rnd = 0d;
        }

        targetSpacing = distanceRange[0] + dx;
        iBeg = 1;
        iEnd = 0;
        xLeast = TauUtil.DMAX;

        // Scan the ray parameter grid looking for points to kill.
        for (int j = 1; j <= i; j++) {
          xCur = calcX(rayParamGrid, tauGrid, distanceRange, j);

          if (sgn * (xCur - targetSpacing) > dx2) {
            // This point looks OK.  See if we have points to kill.
            if (iEnd >= iBeg) {
              for (int k = iBeg; k <= iEnd; k++) {
                keep[k] = false;
              }
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
          for (int k = iBeg; k <= iEnd; k++) {
            keep[k] = false;
          }
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
   * @param rayParamGrid Normalized ray parameter grid
   * @param tauGrid Normalized tau on the same grid
   * @param distanceRange An array of doubles containing the normalized distance at the branch end
   *     points. Note that the length of this array is expected to be two, with the 0 index
   *     containing the minumum range point, and 1 containing the maximum range point.
   * @param i Grid point where the derivative is required
   * @return Distance corresponding to tau(p_i)
   */
  private double calcX(double[] rayParamGrid, double[] tauGrid, double[] distanceRange, int i) {
    if (i == 0) {
      return distanceRange[0];
    } else if (i == rayParamGrid.length - 1) {
      return distanceRange[1];
    } else {
      double h1 = rayParamGrid[i - 1] - rayParamGrid[i];
      double h2 = rayParamGrid[i + 1] - rayParamGrid[i];
      double hh = h1 * h2 * (rayParamGrid[i - 1] - rayParamGrid[i + 1]);
      h1 = Math.pow(h1, 2d);
      h2 = -Math.pow(h2, 2d);

      return -(h2 * tauGrid[i - 1] - (h2 + h1) * tauGrid[i] + h1 * tauGrid[i + 1]) / hh;
    }
  }

  /**
   * Compute the variance of various possible values to keep.
   *
   * @param k0 First trial decimationArray array index
   * @param k1 Second trial decimationArray array index
   * @param k2 Third trial decimationArray array index
   * @param kt Alternate second trial decimationArray array index
   * @return New trial variance of residuals
   */
  private double variance(int k0, int k1, int k2, int kt) {
    double dx1 = Math.abs(decimationArray[k0] - decimationArray[k1]) - targetSpacing;
    double dx2 = Math.abs(decimationArray[k1] - decimationArray[k2]) - targetSpacing;
    double newVar = currentVariance - (Math.pow(dx1, 2d) + Math.pow(dx2, 2d));

    if (kt > k0 && kt < k2) {
      dx1 = Math.abs(decimationArray[k0] - decimationArray[kt]) - targetSpacing;
      dx2 = Math.abs(decimationArray[kt] - decimationArray[k2]) - targetSpacing;
      newVar += Math.pow(dx1, 2d) + Math.pow(dx2, 2d);
      trialNumToKeep = currentNumToKeep;
    } else {
      dx1 = Math.abs(decimationArray[k0] - decimationArray[k2]) - targetSpacing;
      newVar += Math.pow(dx1, 2d);
      trialNumToKeep = currentNumToKeep - 1;
    }

    return newVar;
  }

  /**
   * Compute the variance from scratch for testing purposes.
   *
   * @param keep A boolean array of keep flags. For each element, if true, keep the corresponding
   *     decimationArray value
   * @return A double containing the variance of absolute differences between kept decimationArray
   *     values minus the target difference
   */
  private double doVar(boolean[] keep) {
    int i = 0, currentNumToKeep = 0;
    double currentVariance = 0d;

    for (int j = 1; j < decimationArray.length; j++) {
      if (keep[j]) {
        currentVariance +=
            Math.pow(Math.abs(decimationArray[j] - decimationArray[i]) - targetSpacing, 2d);
        i = j;
        currentNumToKeep++;
      }
    }

    return currentVariance / currentNumToKeep;
  }
}
