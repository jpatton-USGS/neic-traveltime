package gov.usgs.traveltime.tables;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

import gov.usgs.traveltime.ModConvert;
import gov.usgs.traveltime.TauUtil;
import gov.usgs.traveltime.TtStatus;

/**
 * Take care of everything related to one Earth model.
 * 
 * @author Ray Buland
 *
 */
public class EarthModel {
	String earthModel;					// Model name
	ModelSample innerCore;			// Model at the inner core boundary
	ModelSample outerCore;			// Model at the outer core boundary
	ModelSample upperMantle;		// Model at the upper mantle discontinuity
	ModelSample moho;						// Model at the Moho discontinuity
	ModelSample conrad;					// Model at the Conrad discontinuity
	ModelSample surface;				// Model at the free surface
	ArrayList<ModelSample> model;			// Model storage
	ArrayList<ModelShell> shells;			// Model shell parameters
	ModelInterp interp;
	ModConvert convert;
		
	/**
	 * It doesn't take much to get started!
	 * 
	 * @param earthModel Name of the Earth model
	 * @param isCubic True if cubic spline interpolation is to be used
	 */
	public EarthModel(String earthModel, boolean isCubic) {
		this.earthModel = earthModel;
		model = new ArrayList<ModelSample>();
		shells = new ArrayList<ModelShell>();
		interp = new ModelInterp(model, shells, isCubic);
		// Set up the properties.
		if(TauUtil.modelPath == null) {
			TauUtil.getProperties();
		}
	}
	
	/**
	 * Read the Earth model file, set up shells, refine internal boundaries, 
	 * and initialize critical points.
	 * 
	 * @return Travel-time status
	 */
	public TtStatus readModel() {
		String modelFile, modelName;
		boolean disc;
		int n, i = 0, last = 0;
		double rUpperMantle, rMoho, rConrad, rSurface;
		/*
		 * We have to read everything in, but we don't need density, 
		 * anisotropy, or attenuation for the travel-times.
		 */
		@SuppressWarnings("unused")
		double r, rho, vpv, vph, vsv, vsh, eta, qMu, qKappa, rLast = 0d;
		BufferedInputStream inModel;
		Scanner scan;
		
		// Open and read the phase groups file.
		try {
			modelFile = "m"+earthModel+".mod";
			inModel = new BufferedInputStream(new FileInputStream(
					TauUtil.model(modelFile)));
		} catch (FileNotFoundException e) {
			return TtStatus.BAD_MODEL_READ;
		}
		scan = new Scanner(inModel);
		
		// Read the header.
		modelName = scan.next();
		if(!modelName.equals(earthModel)) {
			System.out.println("\n***** Error: model name mismatch (" +
					earthModel + " != " + modelName + ") *****\n");
			scan.close();
			return TtStatus.BAD_MODEL_FILE;
		}
		n = scan.nextInt();
		if(!scan.hasNextInt()) {
			rSurface = scan.nextDouble();
			rUpperMantle = scan.nextDouble();
			rMoho = scan.nextDouble();
			rConrad = scan.nextDouble();
		} else {
			rSurface = TablesUtil.DEFRADIUS;
			rUpperMantle = TablesUtil.DEFUPPERMANTLE;
			rMoho = TablesUtil.DEFMOHO;
			rConrad = TablesUtil.DEFCONRAD;
		}
		
		// Read the model points.
		while(scan.hasNextInt()) {
			i = scan.nextInt();
			if(i != ++last) {
				System.out.format("\n***** Warning: sample %d found, %d " + 
						"expected *****\n\n", last, i);
				last = i;
			}
			r = scan.nextDouble();
			if(r < rLast) {
				System.out.format("\n***** Error: radius %7.2f out of order " + 
						"*****\n\n", r);
				scan.close();
				return TtStatus.BAD_MODEL_FILE;
			}
			rho = scan.nextDouble();
			vpv = scan.nextDouble();
			vph = scan.nextDouble();
			vsv = scan.nextDouble();
			vsh = scan.nextDouble();
			eta = scan.nextDouble();
			qMu = scan.nextDouble();
			qKappa = scan.nextDouble();
			// Trap discontinuities.
			if(r == rLast) {
				disc = true;
				if(model.size() > 0) {
					model.get(model.size()-1).setDisc();
					shells.get(shells.size()-1).addEnd(model.size()-1, r);
				}
				shells.add(new ModelShell(model.size(), r));
			} else {
				disc = false;
			}
			model.add(new ModelSample(r, vpv, vph, vsv, vsh, eta, disc));
			rLast = r;
		}
		// Done, finalize the outermost shell.
		shells.get(shells.size()-1).addEnd(model.size()-1, rLast);
		// Do some crude checks.
		if(i != n) {
			System.out.format("\n***** Warning: %d points found, %d " + 
					"expected *****\n\n", n, i);
		}
		if(rLast != rSurface) {
			System.out.format("\n***** Warning: radius of the model is not the " + 
					"same as the radius of the Earth (%7.2f != %7.2f)\n", rLast, 
					rSurface);
		}
		// OK. We're good (probably).
		scan.close();
		// Set the S velocity to the P velocity in the inner core.
		elimPKJKP();
		// Interpolate velocity.
		interp.interpVel();
		// Find important internal boundaries.
		refineBoundaries(rUpperMantle, rMoho, rConrad, rSurface);
		// Initialize the model specific conversion constants.
		convert = new ModConvert(upperMantle.r, moho.r, conrad.r, surface.r, surface.vs);
		return TtStatus.SUCCESS;
	}
	
	/**
	 * Match the key radii provided with model discontinuities.  This is necessary 
	 * to figure out the phase codes.
	 * 
	 * @param rUpperMantle Preliminary radius of the upper mantle discontinuity 
	 * in kilometers
	 * @param rMoho Preliminary radius of the Moho discontinuity in kilometers
	 * @param rConrad Preliminary radius of the Conrad discontinuity in kilometers
	 * @param rSurface Preliminary radius of the free surface in kilometers
	 */
	private void refineBoundaries(double rUpperMantle, double rMoho, 
			double rConrad, double rSurface) {
		double rInnerCore = TablesUtil.DEFINNERCORE, rOuterCore = TablesUtil.DEFOUTERCORE;
		double tempIC = TauUtil.DMAX, tempOC = TauUtil.DMAX, tempUM = TauUtil.DMAX, 
				tempM = TauUtil.DMAX, tempC = TauUtil.DMAX, tempFS = TauUtil.DMAX;
		double rDisc;
		
		// Find the closest boundary to target boundaries.
		for(int j=0; j<shells.size(); j++) {
			rDisc = shells.get(j).rTop;
			if(Math.abs(rDisc-rInnerCore) < Math.abs(tempIC-rInnerCore)) {
				tempIC = rDisc;
				innerCore = model.get(shells.get(j).iTop);
			}
			if(Math.abs(rDisc-rOuterCore) < Math.abs(tempOC-rOuterCore)) {
				tempOC = rDisc;
				outerCore = model.get(shells.get(j).iTop);
			}
			if(Math.abs(rDisc-rUpperMantle) < Math.abs(tempUM-rUpperMantle)) {
				tempUM = rDisc;
				upperMantle = model.get(shells.get(j).iTop);
			}
			if(Math.abs(rDisc-rMoho) < Math.abs(tempM-rMoho)) {
				tempM = rDisc;
				moho = model.get(shells.get(j).iTop);
			}
			if(Math.abs(rDisc-rConrad) < Math.abs(tempC-rConrad)) {
				tempC = rDisc;
				conrad = model.get(shells.get(j).iTop);
			}
			if(Math.abs(rDisc-rSurface) < Math.abs(tempFS-rSurface)) {
				tempFS = rDisc;
				surface = model.get(shells.get(j).iTop);
			}
		}
		
		// Go around again setting the target ray travel distance in each layer.
		for(int j=0; j<shells.size(); j++) {
			rDisc = shells.get(j).rTop;
			if(rDisc <= rInnerCore) {
				shells.get(j).delX = TablesUtil.DELX[0]/rSurface;
			} else if(rDisc <= rOuterCore) {
				shells.get(j).delX = TablesUtil.DELX[1]/rSurface;
			} else if(rDisc <= rUpperMantle) {
				shells.get(j).delX = TablesUtil.DELX[2]/rSurface;
			} else if(rDisc <= rMoho) {
				shells.get(j).delX = TablesUtil.DELX[3]/rSurface;
			} else if(rDisc <= rConrad) {
				shells.get(j).delX = TablesUtil.DELX[4]/rSurface;
			} else {
				shells.get(j).delX = TablesUtil.DELX[5]/rSurface;
			}
		}
	}
	
	/**
	 * Eliminate the poorly observed PKJKP phase in the inner core by 
	 * replacing the S velocity with the P velocity.
	 */
	private void elimPKJKP() {
		ModelShell shell;
		
		shell = shells.get(0);
		for(int j=shell.iBot; j<=shell.iTop; j++) {
			model.get(j).elimPKJKP();
		}
	}
	
	/**
	 * Getter for ModConvert.
	 * 
	 * @return ModConvert
	 */
	public ModConvert getConvert() {
		return convert;
	}
	
	/**
	 * Interpolate to find V(r).
	 * 
	 * @param type Wave type (P = compressional, S = shear)
	 * @param r Radius in kilometers
	 * @return Velocity in kilometers/second at radius r
	 */
	public double getVel(char type, double r) {
		if(type == 'P') {
			return interp.getVp(r);
		} else {
			return interp.getVs(r);
		}
	}
	
	/**
	 * Interpolate to find V(r) in a particular shell.
	 * 
	 * @param type Wave type (P = compressional, S = shear)
	 * @param shell Shell number
	 * @param r Radius in kilometers
	 * @return Compressional velocity in kilometers/second at radius r
	 */
	public double getVel(char type, int shell, double r) {
		if(type == 'P') {
			return interp.getVp(shell, r);
		} else {
			return interp.getVs(shell, r);
		}
	}
	
	/**
	 * Get the size of the model.
	 * 
	 * @return The number of samples in the Earth model.
	 */
	public int size() {
		return model.size();
	}
	
	/**
	 * Get the jth model depth.
	 * 
	 * @param j Sample index
	 * @return Non-dimensional Earth flattened depth
	 */
	public double getZ(int j) {
		return model.get(j).z;
	}
	
	/**
	 * Get the jth model radius.
	 * 
	 * @param j Sample index
	 * @return Dimensional Earth radius in kilometers
	 */
	public double getR(int j) {
		return model.get(j).r;
	}
	
	/**
	 * Get the jth model slowness.
	 * 
	 * @param type Wave type (P = compressional, S = shear)
	 * @param j Sample index
	 * @return Non-dimensional P-wave slowness
	 */
	public double getSlow(char type, int j) {
		if(type == 'P') {
			return model.get(j).slowP;
		} else {
			return model.get(j).slowS;
		}
	}
	
	/**
	 * Print the model.
	 */
	public void printModel() {
		System.out.format("\n%s %d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n", 
				earthModel, model.size(), innerCore.r, outerCore.r, upperMantle.r, 
				moho.r, conrad.r, surface.r);
		for(int j=0; j<model.size(); j++) {
			model.get(j).printSample(j, false, null);
		}
	}
	
	/**
	 * Print the shell limits.
	 */
	public void printShells() {
		System.out.println("\n\t\tShells:");
		for(int j=0; j<shells.size(); j++) {
			shells.get(j).printShell(j);
		}
	}
}
