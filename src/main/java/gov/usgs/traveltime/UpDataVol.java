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

  int lenNew; // Length of the decimated up-going branch
  double[] pDec; // Decimated up-going branch ray parameters
  double[] tauDec; // Decimated up-going branch tau
  UpDataRef ref;

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
   * Set up volatile copies of data that changes with depth. Note that both P and S models are
   * needed. If this is handling the up-going data for P, the primary model would be for P and the
   * secondary model would be for S.
   *
   * @param ref The up-going reference data source
   * @param primaryModel The primary Earth model data source
   * @param secondaryModel the secondary Earth model data source
   * @param modelConversions Model specific conversions
   */
  public UpDataVol(
      UpDataRef ref,
      ModDataVol primaryModel,
      ModDataVol secondaryModel,
      ModConvert modelConversions) {
    this.ref = ref;
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
    int i;
    double xInt;
    boolean corrTau; // True if upgoingTau needs correcting
    int iUp; // Up-going branch index
    int iBot; // Bottoming depth model index
    double zMax; // Depth of maxSlowness below a low velocity zone

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
    if (-sourceDepth <= TauUtil.DTOL) return;
    // Otherwise, copy the desired data into temporary storage.
    iUp = primaryModel.ref.indexUp[sourceDepthModelIndex];
    //	System.out.println("\t\t\tiUp = "+iUp);
    upgoingRayParams = Arrays.copyOf(ref.pTauUp, ref.upgoingTau[iUp].length);
    upgoingTau = Arrays.copyOf(ref.upgoingTau[iUp], upgoingRayParams.length);
    upgoingDistance = Arrays.copyOf(ref.upgoingDistance[iUp], ref.upgoingDistance[iUp].length);

    // See if we need to correct upgoingTau.
    if (Math.abs(ref.pTauUp[iUp] - maxSlowness) <= TauUtil.DTOL) corrTau = false;
    else corrTau = true;

    maxSlowness = Math.min(maxSlowness, sourceSlowness);
    // Correct the up-going tau values to the exact source depth.
    /*	System.out.println("Partial integrals: "+(float)sourceSlowness+" - "+
    (float)primaryModel.ref.pMod[sourceDepthModelIndex]+"  "+(float)sourceDepth+" - "+
    (float)primaryModel.ref.zMod[sourceDepthModelIndex]); */
    i = 0;
    for (int j = 0; j < upgoingTau.length; j++) {
      if (ref.pTauUp[j] <= maxSlowness) {
        if (corrTau) {
          //		System.out.println("j  p tau (before): "+(j+1)+" "+
          //				(float)ref.pTauUp[j]+" "+(float)upgoingTau[j]);
          upgoingTau[j] -=
              primaryModelTauInt.intLayer(
                  ref.pTauUp[j],
                  sourceSlowness,
                  primaryModel.ref.pMod[sourceDepthModelIndex],
                  sourceDepth,
                  primaryModel.ref.zMod[sourceDepthModelIndex]);
          //		System.out.println("     tau (after): "+(float)upgoingTau[j]+" "+
          //				(float)ref.pXUp[i]);

          // See if we need to correct an end point distance as well.
          if (Math.abs(ref.pTauUp[j] - ref.pXUp[i]) <= TauUtil.DTOL) {
            xInt = primaryModelTauInt.getXLayer();
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
      zMax = primaryModel.findZ(maxSlowness, false);
      iBot = primaryModel.iSource;
      //	System.out.println("\nLVZ integral: "+(float)maxSlowness+" "+sourceDepthModelIndex+" "+
      //			iBot+" "+(float)sourceSlowness+" "+(float)sourceDepth+" "+(float)maxSlowness+
      //			" "+(float)zMax);
      tauIntLVZToSource =
          primaryModelTauInt.intRange(
              maxSlowness,
              sourceDepthModelIndex,
              iBot,
              sourceSlowness,
              sourceDepth,
              maxSlowness,
              zMax);
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
      zMax = secondaryModel.findZ(maxSlowness, true);
      iBot = secondaryModel.iSource;
      //	System.out.println("\nCnv integral: "+(float)maxSlowness+" "+iBot+" "+
      //			(float)maxSlowness+" "+(float)zMax);
      tauIntOtherSurfaceToSource =
          secondaryModelTauInt.intRange(maxSlowness, 0, iBot - 1, maxSlowness, zMax);
      distIntOtherSurfaceToSource = secondaryModelTauInt.getXSum();
      //	System.out.println("tau x = "+(float)tauIntOtherSurfaceToSource+" "+distIntOtherSurfaceToSource);
    } catch (BadDepthException | TauIntegralException e) {
      tauIntOtherSurfaceToSource = 0d;
      distIntOtherSurfaceToSource = 0d;
      //	System.out.println("\nNo Cnv correction needed");
    }
  }

  /**
   * Generate the up-going branch that will be used to compute travel times. The stored up-going
   * branches must be complete in ray parameter samples in order to correct all other travel-time
   * branches to the desired source depth. However, due to the irregular spacing of the ray
   * parameter grid, the interpolation will be unstable. Therefore, the up-going branch must be
   * decimated to be useful later. For very shallow sources, even the decimated grid will be
   * unstable and must be completely replaced.
   *
   * @param pBrn Normalized raw ray parameter grid
   * @param tauBrn Normalized raw tau grid
   * @param xRange Normalized distance range
   * @param xMin Normalized minimum distance interval desired
   * @return A new grid of ray parameter values for the up-going branch
   * @throws TauIntegralException If the tau integration fails
   */
  public double[] realUp(double pBrn[], double tauBrn[], double[] xRange, double xMin)
      throws TauIntegralException {
    int power, len;
    double depth, dp;
    boolean[] keep;

    depth = modelConversions.realZ(sourceDepth);
    if (depth <= modelConversions.zNewUp) {
      // For shallow sources, recompute tau on a more stable ray
      // parameter grid.  The parameters are depth dependent.
      if (depth < 1.5) {
        lenNew = 5;
        power = 6;
      } else if (depth < 10.5) {
        lenNew = 6;
        power = 6;
      } else {
        lenNew = 6;
        power = 7;
      }
      // Allocate some space.
      pDec = new double[lenNew];
      tauDec = new double[lenNew];

      // Create the up-going branch.
      pDec[0] = pBrn[0];
      tauDec[0] = tauBrn[0];
      dp = 0.75d * maxSlowness / Math.pow(lenNew - 2, power);
      for (int j = 1; j < lenNew - 1; j++) {
        pDec[j] = maxSlowness - dp * Math.pow(lenNew - j - 1, power--);
        tauDec[j] =
            primaryModelTauInt.intRange(
                pDec[j], 0, sourceDepthModelIndex - 1, sourceSlowness, sourceDepth);
      }
      pDec[lenNew - 1] = maxSlowness;
      tauDec[lenNew - 1] = tauIntSurfaceToLVZ;
    } else {
      // For deeper sources, it is enough to decimate the ray
      // parameter grid we already have.
      keep = decimator.fastDecimation(pBrn, tauBrn, xRange, xMin);
      if (keep != null) {
        // Do the decimation.
        len = 0;
        for (int k = 0; k < keep.length; k++) {
          if (keep[k]) len++;
        }
        pDec = new double[len];
        tauDec = new double[len];
        for (int k = 0, l = 0; k < keep.length; k++) {
          if (keep[k]) {
            pDec[l] = pBrn[k];
            tauDec[l++] = tauBrn[k];
          }
        }
      } else {
        // We don't need to decimate.
        pDec = Arrays.copyOf(pBrn, pBrn.length);
        tauDec = Arrays.copyOf(tauBrn, tauBrn.length);
      }
    }
    return pDec;
  }

  /**
   * Get the decimated tau values associated with the decimated ray parameter grid.
   *
   * @return Decimated, normalized tau on the decimated ray parameter grid
   */
  public double[] getDecTau() {
    return tauDec;
  }

  /**
   * Print out the up-going branch data corrected for the source depth.
   *
   * @param full If true print the corrected tau array as well.
   */
  public void dumpCorrUp(boolean full) {
    System.out.println("\n     Up-going " + ref.typeUp + " corrected");
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
    System.out.println("\n     Up-going " + ref.typeUp + " decimated");
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
      for (int k = 0; k < tauDec.length; k++) {
        System.out.format("%3d  %8.6f %11.4e\n", k, pDec[k], tauDec[k]);
      }
    }
  }
}
