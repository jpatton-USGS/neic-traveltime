package gov.usgs.traveltime;

import gov.usgs.traveltime.tables.Decimate;
import gov.usgs.traveltime.tables.TauInt;
import gov.usgs.traveltime.tables.TauIntegralException;
import java.util.Arrays;

/**
 * Store volatile up-going branch data for one wave type. Note that all data have been normalized.
 *
 * @author Ray Buland
 */
public class UpDataVol {
  /** An integer containing the source depth model index */
  private int sourceDepthModelIndex;

  /** A double containing the normalized source depth */
  private double sourceDepth;

  /** A double containing the slowness at the source depth */
  private double sourceSlowness;

  /** A double containing the lid slowness if source is in a low velocity zone */
  private double maxSlowness;

  /** An array of doubles containing corrected up-going branch ray parameters */
  private double[] upgoingRayParams;

  /** An array of doubles containing corrected up-going branch tau */
  private double[] upgoingTau;

  /** An array of doubles containing corrected up-going branch distance */
  private double[] upgoingDistance;

  /**
   * A double containing the tau integral from the surface to the low velocity zone for this wave
   * type
   */
  private double tauIntSurfaceToLVZ;

  /**
   * A double containing the tau integral from the low velocity zone to the source for this wave
   * type
   */
  private double tauIntLVZToSource;

  /** A double containing the tau integral from surface to source for other wave type */
  private double tauIntOtherSurfaceToSource;

  /**
   * A double containing the distance integral from the surface to the low velocity zone for this
   * wave type
   */
  private double distIntSurfaceToLVZ;

  /**
   * A double containing the distance integral from the low velocity zone to the source for this
   * wave type
   */
  private double distIntLVZToSource;

  /** A double containing the distance integral from surface to source for other wave type */
  private double distIntOtherSurfaceToSource;

  /** An integer containing the length the decimated up-going branch */
  private int decimatedUpBranchLen;

  /** An array of double values containing the decimated up-going branch ray parameters */
  private double[] decimatedUpBranchRayParams;

  /** An array of double values containing the decimated up-going tau */
  private double[] decimatedUpBranchTau;

  /** An UpDataRef object containing the up data references */
  private UpDataRef upDataReference;

  /** A ModDataVol object containing the primary model */
  private ModDataVol primaryModel;

  /** A ModDataVol object containing the secondary model */
  private ModDataVol secondaryModel;

  /** A ModConvert object containing model specific conversions */
  private ModConvert modelConversions;

  /** A TauInt object containing the primary model tau integration routines */
  private TauInt primaryModelTauInt;

  /** A TauInt object containing the secondary model tau integration routines */
  private TauInt secondaryModelTauInt;

  /** A Decimate object containing the decimation routines */
  Decimate decimator;

  /**
   * Function to return the slowness at the source depth.
   *
   * @return A double containing the slowness at the source depth
   */
  public double getSourceSlowness() {
    return sourceSlowness;
  }

  /**
   * Function to return the lid slowness if source is in a low velocity zone
   *
   * @return A double containing the lid slowness
   */
  public double getMaxSlowness() {
    return maxSlowness;
  }

  /**
   * Function to return the corrected up-going branch ray parameters
   *
   * @return An array of doubles containing corrected up-going branch ray parameters
   */
  public double[] getUpgoingRayParams() {
    return upgoingRayParams;
  }

  /**
   * Function to return the corrected up-going branch tau
   *
   * @return An array of doubles containing corrected up-going branch tau
   */
  public double[] getUpgoingTau() {
    return upgoingTau;
  }

  /**
   * Function to return the corrected up-going branch distance
   *
   * @return An array of doubles containing corrected up-going branch distance
   */
  public double[] getUpgoingDistance() {
    return upgoingDistance;
  }

  /**
   * Function to return the tau integral from the surface to the low velocity zone
   *
   * @return A double containing the tau integral from the surface to the low velocity zone for this
   *     wave type
   */
  public double getTauIntSurfaceToLVZ() {
    return tauIntSurfaceToLVZ;
  }

  /**
   * Function to return the tau integral from the low velocity zone to the source
   *
   * @return A double containing the tau integral from tlow velocity zone to the source for this
   *     wave type
   */
  public double getTauIntLVZToSource() {
    return tauIntLVZToSource;
  }

  /**
   * Function to return the tau integral from surface to source for other wave type
   *
   * @return A double containing the tau integral from surface to source for other wave type
   */
  public double getTauIntOtherSurfaceToSource() {
    return tauIntOtherSurfaceToSource;
  }

  /**
   * Function to return the distance integral from the surface to the low velocity zone
   *
   * @return A double containing the distance integral from the surface to the low velocity zone for
   *     this wave type
   */
  public double getDistIntSurfaceToLVZ() {
    return distIntSurfaceToLVZ;
  }

  /**
   * Function to return the distance integral from the low velocity zone to the source
   *
   * @return A double containing the distance integral from tlow velocity zone to the source for
   *     this wave type
   */
  public double getDistIntLVZToSource() {
    return distIntLVZToSource;
  }

  /**
   * Function to return the distance integral from surface to source for other wave type
   *
   * @return A double containing the distance integral from surface to source for other wave type
   */
  public double getDistIntOtherSurfaceToSource() {
    return distIntOtherSurfaceToSource;
  }

  /**
   * Function to return the up data references
   *
   * @return An UpDataRef object containing the up data references
   */
  public UpDataRef getUpDataReference() {
    return upDataReference;
  }

  /**
   * Function to return the tau values associated with the decimated ray parameter grid.
   *
   * @return An array of doubles containing the decimated, normalized tau on the decimated ray
   *     parameter grid
   */
  public double[] getDecimatedUpBranchTau() {
    return decimatedUpBranchTau;
  }

  /**
   * Set up volatile copies of data that changes with depth. Note that both P and S models are
   * needed. If this is handling the up-going data for P, the primary model would be for P and the
   * secondary model would be for S.
   *
   * @param upDataReference The up-going reference data source
   * @param primaryModel The primary Earth model data source
   * @param secondaryModel the secondary Earth model data source
   * @param modelConversions Model specific conversions
   */
  public UpDataVol(
      UpDataRef upDataReference,
      ModDataVol primaryModel,
      ModDataVol secondaryModel,
      ModConvert modelConversions) {
    this.upDataReference = upDataReference;
    this.primaryModel = primaryModel;
    this.secondaryModel = secondaryModel;
    this.modelConversions = modelConversions;

    // Set up the integration routines.
    primaryModelTauInt = new TauInt(primaryModel);
    secondaryModelTauInt = new TauInt(secondaryModel);

    // Set up the decimation.
    decimator = new Decimate();
  }

  /**
   * Correct up-going tau to the desired source depth. The up-going branches are used to correct tau
   * for all ray parameters for all travel-time branches. At the same time, integrals are computed
   * for the largest ray parameter (usually equal to the source depth slowness) needed to correct
   * tau for the largest ray parameter for all branches.
   *
   * @param depth Normalized source depth
   * @throws BadDepthException If the source depth is too deep
   * @throws TauIntegralException If the tau integral fails
   */
  public void newDepth(double depth) throws BadDepthException, TauIntegralException {
    int depthModelBottomIndex; // Bottoming depth model index
    double maxSlownessDepth; // Depth of maxSlowness below a low velocity zone

    // Initialize.
    sourceDepth = depth;
    tauIntSurfaceToLVZ = 0d;
    tauIntLVZToSource = 0d;
    tauIntOtherSurfaceToSource = 0d;
    distIntSurfaceToLVZ = 0d;
    distIntLVZToSource = 0d;
    distIntOtherSurfaceToSource = 0d;

    // Get the source slowness.
    sourceSlowness = primaryModel.findP(sourceDepth);
    sourceDepthModelIndex = primaryModel.iSource;
    maxSlowness = primaryModel.findMaxP();
    //	primaryModel.printFind(false);

    // If the source is at the surface, we're already done.
    if (-sourceDepth <= TauUtil.DTOL) {
      return;
    }

    // Otherwise, copy the desired data into temporary storage.
    int upgoingBranchIndex = primaryModel.ref.indexUp[sourceDepthModelIndex];
    //	System.out.println("\t\t\tiUp = "+upgoingBranchIndex);
    upgoingRayParams =
        Arrays.copyOf(
            upDataReference.pTauUp, upDataReference.upgoingTau[upgoingBranchIndex].length);
    upgoingTau =
        Arrays.copyOf(upDataReference.upgoingTau[upgoingBranchIndex], upgoingRayParams.length);
    upgoingDistance =
        Arrays.copyOf(
            upDataReference.upgoingDistance[upgoingBranchIndex],
            upDataReference.upgoingDistance[upgoingBranchIndex].length);

    // See if we need to correct upgoingTau
    boolean correctTau = true; // True if upgoingTau needs correcting
    if (Math.abs(upDataReference.pTauUp[upgoingBranchIndex] - maxSlowness) <= TauUtil.DTOL) {
      correctTau = false;
    } else {
      correctTau = true;
    }

    maxSlowness = Math.min(maxSlowness, sourceSlowness);
    // Correct the up-going tau values to the exact source depth.
    /*	System.out.println("Partial integrals: "+(float)sourceSlowness+" - "+
    (float)primaryModel.ref.pMod[sourceDepthModelIndex]+"  "+(float)sourceDepth+" - "+
    (float)primaryModel.ref.zMod[sourceDepthModelIndex]); */
    int i = 0;
    for (int j = 0; j < upgoingTau.length; j++) {
      if (upDataReference.pTauUp[j] <= maxSlowness) {
        if (correctTau) {
          //		System.out.println("j  p tau (before): "+(j+1)+" "+
          //				(float)upDataReference.pTauUp[j]+" "+(float)upgoingTau[j]);
          upgoingTau[j] -=
              primaryModelTauInt.intLayer(
                  upDataReference.pTauUp[j],
                  sourceSlowness,
                  primaryModel.ref.pMod[sourceDepthModelIndex],
                  sourceDepth,
                  primaryModel.ref.zMod[sourceDepthModelIndex]);
          //		System.out.println("     tau (after): "+(float)upgoingTau[j]+" "+
          //				(float)upDataReference.pXUp[i]);

          // See if we need to correct an end point distance as well.
          if (Math.abs(upDataReference.pTauUp[j] - upDataReference.pXUp[i]) <= TauUtil.DTOL) {
            double xInt = primaryModelTauInt.getXLayer();
            upgoingDistance[i++] -= xInt;
            //			System.out.println("i  x (after) dx = "+i+" "+
            //					(float)upgoingDistance[i-1]+" "+(float)xInt);
          }
        }
      } else break;
    }

    /**
     * Compute tau and distance for the ray parameter equal to the source slowness (i.e., horizontal
     * take-off angle from the source).
     */
    //	System.out.println("\nEnd integral: "+(float)maxSlowness+" "+sourceDepthModelIndex+" "+
    //			(float)sourceSlowness+" "+(float)sourceDepth);
    tauIntSurfaceToLVZ =
        primaryModelTauInt.intRange(
            maxSlowness, 0, sourceDepthModelIndex - 1, sourceSlowness, sourceDepth);
    distIntSurfaceToLVZ = primaryModelTauInt.getXSum();
    //	System.out.println("tau x = "+(float)tauIntSurfaceToLVZ+" "+distIntSurfaceToLVZ);

    /**
     * If the source depth is in a low velocity zone, we need to compute tau and distance down to
     * the shallowest turning ray (the horizontal ray is trapped).
     */
    if (maxSlowness > sourceSlowness) {
      maxSlownessDepth = primaryModel.findZ(maxSlowness, false);
      depthModelBottomIndex = primaryModel.iSource;
      //	System.out.println("\nLVZ integral: "+(float)maxSlowness+" "+sourceDepthModelIndex+" "+
      //			depthModelBottomIndex+" "+(float)sourceSlowness+" "+(float)sourceDepth+"
      // "+(float)maxSlowness+
      //			" "+(float)maxSlownessDepth);
      tauIntLVZToSource =
          primaryModelTauInt.intRange(
              maxSlowness,
              sourceDepthModelIndex,
              depthModelBottomIndex,
              sourceSlowness,
              sourceDepth,
              maxSlowness,
              maxSlownessDepth);
      distIntLVZToSource = primaryModelTauInt.getXSum();
      //	System.out.println("tau x = "+(float)tauIntLVZToSource+" "+distIntLVZToSource);
    } else {
      tauIntLVZToSource = 0d;
      distIntLVZToSource = 0d;
    }

    /**
     * Compute tau and distance for the other wave type for the ray parameter equal to the source
     * slowness.
     */
    try {
      maxSlownessDepth = secondaryModel.findZ(maxSlowness, true);
      depthModelBottomIndex = secondaryModel.iSource;
      //	System.out.println("\nCnv integral: "+(float)maxSlowness+" "+depthModelBottomIndex+" "+
      //			(float)maxSlowness+" "+(float)maxSlownessDepth);
      tauIntOtherSurfaceToSource =
          secondaryModelTauInt.intRange(
              maxSlowness, 0, depthModelBottomIndex - 1, maxSlowness, maxSlownessDepth);
      distIntOtherSurfaceToSource = secondaryModelTauInt.getXSum();
      //	System.out.println("tau x = "+(float)tauIntOtherSurfaceToSource+"
      // "+distIntOtherSurfaceToSource);
    } catch (BadDepthException | TauIntegralException e) {
      tauIntOtherSurfaceToSource = 0d;
      distIntOtherSurfaceToSource = 0d;
      //	System.out.println("\nNo Cnv correction needed");
    }
  }

  /**
   * Function to generate the up-going branch that will be used to compute travel times. The stored
   * up-going branches must be complete in ray parameter samples in order to correct all other
   * travel-time branches to the desired source depth. However, due to the irregular spacing of the
   * ray parameter grid, the interpolation will be unstable. Therefore, the up-going branch must be
   * decimated to be useful later. For very shallow sources, even the decimated grid will be
   * unstable and must be completely replaced.
   *
   * @param normalizedRawRayParams An array of doubles containing the normalized raw ray parameter
   *     grid
   * @param normalizedRawTau An array of doubles containing the normalized raw tau grid
   * @param normalizedDistanceRange An array of doubles containing the normalized distance range
   * @param normalizedMinDistInterval A double value containing the normalized minimum distance
   *     interval desired
   * @return An array of doubles containing the new grid of ray parameter values for the up-going
   *     branch
   * @throws TauIntegralException If the tau integration fails
   */
  public double[] realUp(
      double normalizedRawRayParams[],
      double normalizedRawTau[],
      double[] normalizedDistanceRange,
      double normalizedMinDistInterval)
      throws TauIntegralException {
    int power, len;
    double depth, dp;
    boolean[] keep;

    depth = modelConversions.realZ(sourceDepth);
    if (depth <= modelConversions.zNewUp) {
      // For shallow sources, recompute tau on a more stable ray
      // parameter grid.  The parameters are depth dependent.
      if (depth < 1.5) {
        decimatedUpBranchLen = 5;
        power = 6;
      } else if (depth < 10.5) {
        decimatedUpBranchLen = 6;
        power = 6;
      } else {
        decimatedUpBranchLen = 6;
        power = 7;
      }

      // Allocate some space.
      decimatedUpBranchRayParams = new double[decimatedUpBranchLen];
      decimatedUpBranchTau = new double[decimatedUpBranchLen];

      // Create the up-going branch.
      decimatedUpBranchRayParams[0] = normalizedRawRayParams[0];
      decimatedUpBranchTau[0] = normalizedRawTau[0];
      dp = 0.75d * maxSlowness / Math.pow(decimatedUpBranchLen - 2, power);

      for (int j = 1; j < decimatedUpBranchLen - 1; j++) {
        decimatedUpBranchRayParams[j] =
            maxSlowness - dp * Math.pow(decimatedUpBranchLen - j - 1, power--);
        decimatedUpBranchTau[j] =
            primaryModelTauInt.intRange(
                decimatedUpBranchRayParams[j],
                0,
                sourceDepthModelIndex - 1,
                sourceSlowness,
                sourceDepth);
      }
      decimatedUpBranchRayParams[decimatedUpBranchLen - 1] = maxSlowness;
      decimatedUpBranchTau[decimatedUpBranchLen - 1] = tauIntSurfaceToLVZ;
    } else {
      // For deeper sources, it is enough to decimate the ray
      // parameter grid we already have.
      keep =
          decimator.fastDecimation(
              normalizedRawRayParams,
              normalizedRawTau,
              normalizedDistanceRange,
              normalizedMinDistInterval);

      if (keep != null) {
        // Do the decimation.
        len = 0;

        for (int k = 0; k < keep.length; k++) {
          if (keep[k]) {
            len++;
          }
        }

        decimatedUpBranchRayParams = new double[len];
        decimatedUpBranchTau = new double[len];

        for (int k = 0, l = 0; k < keep.length; k++) {
          if (keep[k]) {
            decimatedUpBranchRayParams[l] = normalizedRawRayParams[k];
            decimatedUpBranchTau[l++] = normalizedRawTau[k];
          }
        }
      } else {
        // We don't need to decimate.
        decimatedUpBranchRayParams =
            Arrays.copyOf(normalizedRawRayParams, normalizedRawRayParams.length);
        decimatedUpBranchTau = Arrays.copyOf(normalizedRawTau, normalizedRawTau.length);
      }
    }
    return decimatedUpBranchRayParams;
  }

  /**
   * Print out the up-going branch data corrected for the source depth.
   *
   * @param full If true print the corrected tau array as well.
   */
  public void dumpCorrUp(boolean full) {
    System.out.println("\n     Up-going " + upDataReference.typeUp + " corrected");
    System.out.format(
        "TauEnd: %8.6f %8.6f %8.6f  XEnd: %8.6f %8.6f %8.6f\n",
        tauIntSurfaceToLVZ,
        tauIntLVZToSource,
        tauIntOtherSurfaceToSource,
        distIntSurfaceToLVZ,
        distIntLVZToSource,
        distIntOtherSurfaceToSource);
    if (full) {
      System.out.println("          p        tau");
      for (int k = 0; k < upgoingTau.length; k++) {
        System.out.format("%3d  %8.6f %11.4e\n", k, upgoingRayParams[k], upgoingTau[k]);
      }
      /*	if(brnLen > upgoingTau.length) {
      	System.out.format("%3d  %8.6f  %8.6f\n",brnLen-1,upgoingRayParams[brnLen-1],
      			tauIntSurfaceToLVZ+tauIntLVZToSource);
      } */
    }
  }

  /**
   * Print out the decimated up-going branch data corrected for the source depth.
   *
   * @param full If true print the corrected tau array as well.
   */
  public void dumpDecUp(boolean full) {
    System.out.println("\n     Up-going " + upDataReference.typeUp + " decimated");
    System.out.format(
        "TauEnd: %8.6f %8.6f %8.6f  XEnd: %8.6f %8.6f %8.6f\n",
        tauIntSurfaceToLVZ,
        tauIntLVZToSource,
        tauIntOtherSurfaceToSource,
        distIntSurfaceToLVZ,
        distIntLVZToSource,
        distIntOtherSurfaceToSource);
    if (full) {
      System.out.println("          p        tau");
      for (int k = 0; k < decimatedUpBranchTau.length; k++) {
        System.out.format(
            "%3d  %8.6f %11.4e\n", k, decimatedUpBranchRayParams[k], decimatedUpBranchTau[k]);
      }
    }
  }
}
