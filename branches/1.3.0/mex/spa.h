/////////////////////////////////////////////
// HEADER FILE for SPA.C //
// //
// Solar Position Algorithm (SPA) //
// for //
// Solar Radiation Application //
// //
// May 12, 2003 //
// //
// Filename: SPA.H //
// //
// Afshin Michael Andreas //
// afshin_andreas@nrel.gov (303)384-6383 //
// //
// Measurement & Instrumentation Team //
// Solar Radiation Research Laboratory //
// National Renewable Energy Laboratory //
// 1617 Cole Blvd, Golden, CO 80401 //
/////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// //
// Usage: //
// //
// 1) In calling program, include this header file, //
// by adding this line to the top of file: //
// #include "spa.h" //
// //
// 2) In calling program, declare the SPA structure: //
// spa_data spa; //
// //
// 3) Enter the required input values into SPA structure //
// (input values listed in comments below) //
// //
// 4) Call the SPA calculate function and pass the SPA structure //
// (prototype is declared at the end of this header file): //
// spa_calculate(&spa); //
// //
// All output values (listed in comments below) will be //
// computed and returned in the passed SPA structure. //
// //
////////////////////////////////////////////////////////////////////////
#ifndef __solar_position_algorithm_header
#define __solar_position_algorithm_header
typedef struct {
	//----------------------INPUT VALUES------------------------
	int year; //4-digit year
	int month; //2-digit month (1-12)
	int day; //2-digit day (1-31)
	int hour; //observer local hour
	int minute; //observer local minute
	int second; //observer local second
	float delta_t; //difference between earth rotation time and terrestrial time
			//(from observation) [seconds]
	float timezone; //observer timezone (negative west of greenwich) [hours]
	float longitude; //observer longitude (negative west of greenwich) [degrees]
	float latitude; //observer latitude (negative south of equator) [degrees]
	float elevation; //observer elevation [meters]
	float pressure; //annual average local pressure [millibars]
	float temperature; //annual average local temperature [degrees celcius]
	float slope; //surface slope (measured from the horizontal plane) [degrees]
	float azm_rotation; //surface azimuth rotation (measured from south to projection of
			//surface normal on horizontal plane, negative west) [degrees]
	float atmos_refract; //atmospheric refraction at sunrise and sunset [degrees]
			//0.5667 degrees is typical
	//-----------------Intermediate OUTPUT VALUES--------------------
	double jd; //julian day
	double jc; //julian century;
	double jde; //julian ephemeris day
	double jce; //julian ephemeris century
	double jme; //julian ephemeris millennium
	double l; //earth heliocentric longitude [degrees]
	double b; //earth heliocentric latitude [degrees]
	double r; //earth radius vector [Astronomical Units, AU]
	double theta; //geocentric longitude [degrees]
	double beta; //geocentric latitude [degrees]
	double x0; //mean elongation (moon-sun) [degrees]
	double x1; //mean anomaly (sun) [degrees]
	double x2; //mean anomaly (moon) [degrees]
	double x3; //argument latitude (moon) [degrees]
	double x4; //ascending longitude (moon) [degrees]
	double del_psi; //nutation longitude [degrees]
	double del_epsilon; //nutation obliquity [degrees]
	double epsilon0; //ecliptic mean obliquity [arc seconds]
	double epsilon; //ecliptic true obliquity [degrees]
	double del_tau; //aberration correction [degrees]
	double lamda; //apparent sun longitude [degrees]
	double nu0; //greenwich mean sidereal time [degrees]
	double nu; //greenwich sidereal time [degrees]
	double alpha; //geocentric sun right ascension [degrees]
	double delta; //geocentric sun declination [degrees]
	double h; //observer hour angle [degrees]
	double xi; //sun equatorial horizontal parallax [degrees]
	double del_alpha; //sun right ascension parallax [degrees]
	double delta_prime; //topocentric sun declination [degrees]
	double alpha_prime; //topocentric sun right ascension [degrees]
	double h_prime; //topocentric local hour angle [degrees]
	double e0; //topocentric elevation angle (uncorrected) [degrees]
	double del_e; //atmospheric refraction correction [degrees]
	double e; //topocentric elevation angle (corrected) [degrees]
	//---------------------Final OUTPUT VALUES------------------------
	double zenith; //topocentric zenith angle [degrees]
	double azimuth180; //topocentric azimuth angle (westward from south) [-180 to 180 degrees]
	double azimuth; //topocentric azimuth angle (eastward from north) [0 to 360 degrees]
	double incidence; //surface incidence angle [degrees]
	double eot; //equation of time [degrees]
	double suntransit; //local sun transit time (or solar noon) [fractional hour]
	double sunrise; //local sunrise time [fractional hour]
	double sunset; //local sunset time [fracitonal hour]
} spa_data;
//Calculate all SPA output values (in structure) based on input values passed in structure
void spa_calculate(spa_data *spa);
#endif
