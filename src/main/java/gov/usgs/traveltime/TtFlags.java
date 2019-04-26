package gov.usgs.traveltime;

/**
 * Organize all the auxiliary travel-time information by phase.
 * 
 * @author Ray Buland
 *
 */
public class TtFlags {
	final String phGroup;				// Phase group
	final String auxGroup;			// Auxiliary phase group
	final boolean isRegional;		// True if phase is regional
	final boolean isDepth;			// True if phase is depth sensitive
	final boolean canUse;				// True if phase can be used in a location
	final boolean dis;					// True if phase is to be down weighted
	final TtStat ttStat;				// Phase statistics
	final Ellip ellip;					// Ellipticity correction
	final Ellip upEllip;				// Ellipticity correction for up-going P and S
	
	/**
	 * Each phase is created all at once.
	 * 
	 * @param phGroup Phase group
	 * @param auxGroup Auxiliary (complementary) phase group
	 * @param isRegional True if the phase is flagged as regional
	 * @param isDepth True if the phase is depth sensitive
	 * @param canUse True if the phase may be used in an earthquake 
	 * location
	 * @param dis True if the phase should be down weighted
	 * @param ttStat Travel-time statistics object
	 * @param ellip Ellipticity correction object
	 * @param upEllip Ellipticity correction object for an associated 
	 * up going branch if any
	 */
	public TtFlags(String phGroup, String auxGroup, boolean isRegional, 
			boolean isDepth, boolean canUse, boolean dis, TtStat ttStat, 
			Ellip ellip, Ellip upEllip) {
		this.phGroup = phGroup;
		this.auxGroup = auxGroup;
		this.isRegional = isRegional;
		this.isDepth = isDepth;
		this.canUse = canUse;
		this.dis = dis;
		this.ttStat = ttStat;
		this.ellip = ellip;
		this.upEllip = upEllip;
	}
}