function [Head,Data]=gFXtfRead000(XtfHead,SubCh)
%Read [Head,Data] from XtfHead.fName (*.xtf) file for Message Type 000 (Sonar Data Message).
%
% XsfHead - Xsf Header structure;
% Head - Header structure;
% Data - Data Body for sonar channels.
% Head include the next addition fields: Head.HSubChannelNumber, Head.HMessageNum.
% Example: [Head,Data]=gFXtfRead000(XtfHead,0);
%
% mail@ge0mlib.ru 01/08/2016
%
% This function comes from geomlib (FEX: 58387, BSD Licence)

% $Id: gFXtfRead000.m 11303 2018-05-28 21:39:31Z Joaquim Luis $

	[fId, mes]=fopen(XtfHead.fName,'r');
	if ~isempty(mes), error(['gFXtfRead000: ' mes]),	end
	LHead = (XtfHead.RHeaderType==0) & (XtfHead.RSubChannelNumber==SubCh);
	LenHead=sum(LHead);		nHead=find(LHead);
	ChF=max(XtfHead.RNumChansToFollow(nHead));

	% Begin Header and ChanInfo Allocate
	v(1,LenHead) = single(0);
	vc(ChF,LenHead) = single(0);
	Head=struct('HSubChannelNumber',SubCh,'HMessageNum',v,'HYear',v,'HMonth',v,'HDay',v,'HHour',v,'HMinute',v,...
		'HSecond',v,'HHSeconds',v,'HJulianDay',v,'HEventNumber',v,'HPingNumber',v,'HSoundVelocity',v,...
		'HOceanTide',v,'HWaterTemperature',v,'HComputedSoundVelocity',v,'HShipSpeed',v,...
		'HShipGyro',v,'HShipYcoordinate',v,'HShipXcoordinate',v,'HShipAltitude',v,'HShipDepth',v,...
		'HFixTimeHour',v,'HFixTimeMinute',v,'HFixTimeSecond',v,'HFixTimeHsecond',v,'HSensorSpeed',v,...
		'HSensorYcoordinate',v,'HSensorXcoordinate',v,'HSonarStatus',v,'HRangeToFish',v,'HBearingToFish',v,...
		'HSensorDepth',v,'HSensorPrimaryAltitude',v,...
		'HSensorAuxAltitude',v,'HSensorPitch',v,'HSensorRoll',v,'HSensorHeading',v,'HHeave',v,'HYaw',v,...
		'CChannelNumber',vc,'CSlantRange',vc,'CGroundRange',vc,'CNumSamples',vc);

	df=0;	fseek(fId,0,'bof');	CDataSeek=vc;
	for (n = 1:LenHead)
		fseek(fId,XtfHead.RSeek(nHead(n))-df,'cof');
		% Begin Header Read
		Head.HMessageNum(n)=nHead(n);
		Head.HYear(n)=fread(fId,1,'uint16');					%Ping year
		Head.HMonth(n)=fread(fId,1,'uint8');					%Ping month
		Head.HDay(n)=fread(fId,1,'uint8');						%Ping day
		Head.HHour(n)=fread(fId,1,'uint8');						%Ping hour
		Head.HMinute(n)=fread(fId,1,'uint8');					%Ping minute
		Head.HSecond(n)=fread(fId,1,'uint8');					%Ping seconds
		Head.HHSeconds(n)=fread(fId,1,'uint8');					%Ping hundredths of seconds (0-99)
		Head.HJulianDay(n)=fread(fId,1,'uint16');				%Julian day of a ping’s occurrence.
		Head.HEventNumber(n)=fread(fId,1,'uint32');				%Last logged event number; nav interface template token=O. NOTE: In Isis v4.30 and earlier EventNumber field was located at byte 26-27 and was a two byte WORD.  At byte 24-25 there used to be a WORD CurrentLineID.  The CurrentLineID field no longer exists in the .XTF format.  Therefore, to read the event number correctly an application MUST check the Isis version string starting at byte 10 of the XTFFILEHEADER structure.
		Head.HPingNumber(n)=fread(fId,1,'uint32');				%Counts consecutively (usually from 0) and increments for each update.  Isis Note: The counters are different between sonar and bathymetry updates.
		Head.HSoundVelocity(n)=fread(fId,1,'float32');			%m/s,  Isis uses 750 (one way), some XTF files use 1500.  Note: Can be changed on Isis menu.  This value is never computed and can only be changed manually by the user. See ComputedSoundVelocity below.
		Head.HOceanTide(n)=fread(fId,1,'float32');				%Altitude above Geoide (from RTK), if present; ELSE Ocean tide in meters; nav interface template token = {t} Isis Note: Can be changed by the user on the Configure menu in Isis.
 		%fseek(fId, 4, 0);	%Head.HReserved2(n)=fread(fId,1,'uint32');		%Unused. Set to 0.
		%fseek(fId, 4, 0);	%Head.HConductivityFreq(n)=fread(fId,1,'float32');		%Conductivity frequency in Hz. nav interface template token = Q Raw CTD information.  The Freq values are those sent up by the Seabird CTD. The Falmouth Scientific CTD sends up computed data.
		%fseek(fId, 4, 0);	%Head.HTemperatureFreq(n)=fread(fId,1,'float32');		%Temperature frequency in Hz. nav interface template token = b Raw CTD information.  The Freq values are those sent up by the Seabird CTD. The Falmouth Scientific CTD sends up computed data.
		%fseek(fId, 4, 0);	%Head.HPressureFreq(n)=fread(fId,1,'float32');			%Pressure frequency in Hz. nav interface template token = 0.  Raw CTD information.  The Freq values are those sent up by the Seabird CTD. The Falmouth Scientific CTD sends up computed data.
		%fseek(fId, 4, 0);	%Head.HPressureTemp(n)=fread(fId,1,'float32');			%Pressure temperature (Degrees C); nav interface template token = ; Raw CTD information.  The Freq values are those sent up by the Seabird CTD. The Falmouth Scientific CTD sends up computed data.
		%fseek(fId, 4, 0);	%Head.HConductivity(n)=fread(fId,1,'float32');			%Conductivity in Siemens/m; nav interface token = {c}; can be computed from Q Computed CTD information. When using a Seabird CTD, these values are computed from the raw Freq values (above).
		fseek(fId, 24, 0);
		Head.HWaterTemperature(n)=fread(fId,1,'float32');		%Water temperature in Celsius. nav interface token = {w}; can be computed from b. Computed CTD information. When using a Seabird CTD, these values are computed from the raw Freq values (above)
		fseek(fId, 4, 0);	%Head.HPressure(n)=fread(fId,1,'float32');				%Water pressure in psia; nav interface token = {p}; can be computed from 0. Computed CTD information. When using a Seabird CTD, these values are computed from the raw Freq values (above).
		Head.HComputedSoundVelocity(n)=fread(fId,1,'float32');	%Meters/second computed from Conductivity, WaterTemperature, and Pressure using the Chen Millero formula (1977), formula (JASA, 62, 1129-1135)
		%fseek(fId, 4, 0);	%Head.HMagX(n)=fread(fId,1,'float32');			%X-axis magnetometer data in mgauss. Nav interface template token = e. Sensors Information.
		%fseek(fId, 4, 0);	%Head.HMagY(n)=fread(fId,1,'float32');			%Y-axis magnetometer data in mgauss. Nav interface template token = w. Sensors Information.
		%fseek(fId, 4, 0);	%Head.HMagZ(n)=fread(fId,1,'float32');			%Z-axis magnetometer data in mgauss. Nav interface template token = z. Sensors Information.
		%fseek(fId, 4, 0);	%Head.HAuxVal1(n)=fread(fId,1,'float32');		%Sensors Information. Nav interface template token = 1. Auxiliary values can be used to store and display any value at the user's discretion. Not used in any calculation in Isis or Target. Isis Note: Displayed in the “Sensors” window by selecting “Window?Text?Sensors”
		%fseek(fId, 4, 0);	%Head.HAuxVal2(n)=fread(fId,1,'float32');		%Sensors Information. Nav interface template token = 2. Auxiliary values can be used to store and display any value at the user's discretion.   These are not used in any calculation in Isis or Target. Isis Note: Displayed in the “Sensors” window by selecting “Window?Text?Sensors”
		%fseek(fId, 4, 0);	%Head.HAuxVal3(n)=fread(fId,1,'float32');		%Sensors Information. Nav interface template token = 3. Auxiliary values can be used to store and display any value at the user's discretion.   These are not used in any calculation in Isis or Target. Isis Note: Displayed in the “Sensors” window by selecting “Window?Text?Sensors”
		%fseek(fId, 4, 0);	%Head.HAuxVal4(n)=fread(fId,1,'float32');		%Sensors Information. Nav interface template token = 4. Auxiliary values can be used to store and display any value at the user's discretion.   These are not used in any calculation in Isis or Target. Isis Note: Displayed in the “Sensors” window by selecting “Window?Text?Sensors”
		%fseek(fId, 4, 0);	%Head.HAuxVal5(n)=fread(fId,1,'float32');		%Sensors Information. Nav interface template token = 5. Auxiliary values can be used to store and display any value at the user's discretion.   These are not used in any calculation in Isis or Target. Isis Note: Displayed in the “Sensors” window by selecting “Window?Text?Sensors”
		%fseek(fId, 4, 0);	%Head.HAuxVal6(n)=fread(fId,1,'float32');		%Sensors Information. Nav interface template token = 6. Auxiliary values can be used to store and display any value at the user's discretion.   These are not used in any calculation in Isis or Target. Isis Note: Displayed in the “Sensors” window by selecting “Window?Text?Sensors”
		%fseek(fId, 4, 0);	%Head.HSpeedLog(n)=fread(fId,1,'float32');				%Sensors Information. Speed log sensor on towfish in knots; Note: This is not fish speed. Nav interface template token = s.
		%fseek(fId, 4, 0);	%Head.HTurbidity(n)=fread(fId,1,'float32');		%Sensors Information. Turbidity sensor (0 to +5 volts) multiplied by 10000. nav interface template token = | (the “pipe” symbol).
		fseek(fId, 44, 0);
		Head.HShipSpeed(n)=fread(fId,1,'float32');				%Ship Navigation information. Ship speed in knots. nav interface template token = v. Isis Note: These values are stored only and are not part of any equation or computation in Isis.
		Head.HShipGyro(n)=fread(fId,1,'float32');				%Ship Navigation information. Ship gyro in degrees. nav interface template token = G. Isis Note: This is used as the directional sensor for Multibeam Bathymetry data.
		Head.HShipYcoordinate(n)=fread(fId,1,'float64');		%Ship Navigation information. Ship latitude or northing in degrees. nav interface template token = y. Isis Note: These values are stored only and are not part of any equation or computation in Isis.
		Head.HShipXcoordinate(n)=fread(fId,1,'float64');		%Ship Navigation information. Ship longitude or easting in degrees. nav interface template token = x. Isis Note: These values are stored only and are not part of any equation or computation in Isis.
		Head.HShipAltitude(n)=fread(fId,1,'uint16');			%Ship altitude in decimeters
		Head.HShipDepth(n)=fread(fId,1,'uint16');				%Ship depth in decimeters.
		Head.HFixTimeHour(n)=fread(fId,1,'uint8');				%Sensor Navigation information. Hour of most recent nav update. nav interface template token = H. Isis Note: The time of the nav is adjusted by the NavLatency stored in the XTF file header.
		Head.HFixTimeMinute(n)=fread(fId,1,'uint8');			%Sensor Navigation information. Minute of most recent nav update. nav interface template token = I. Isis Note: The time of the nav is adjusted by the NavLatency stored in the XTF file header.
		Head.HFixTimeSecond(n)=fread(fId,1,'uint8');			%Sensor Navigation information. Second of most recent nav update. nav interface template token = S. Isis Note: The time of the nav is adjusted by the NavLatency stored in the XTF file header.
		Head.HFixTimeHsecond(n)=fread(fId,1,'uint8');			%Sensor Navigation information. Hundredth of a Second of most recent nav update. Isis Note: The time of the nav is adjusted by the NavLatency stored in the XTF file header.
		Head.HSensorSpeed(n)=fread(fId,1,'float32');			%Sensor Navigation information. Speed of towfish in knots. Used for speed correction and position calculation; nav interface template token = V.
		fseek(fId, 4, 0);	%Head.HKP(n)=fread(fId,1,'float32');		%Sensor Navigation information. Kilometers Pipe; nav interface template token = {K}.
		Head.HSensorYcoordinate(n)=fread(fId,1,'float64');		%Sensor Navigation information. Sensor latitude or northing; nav interface template token = E. Note: when NavUnits in the file header is 0, values are in meters (northings and eastings).  When NavUnits is 3, values are in Lat/Long.  Also see the Layback value, below.
		Head.HSensorXcoordinate(n)=fread(fId,1,'float64');		%Sensor Navigation information. Sensor longitude or easting; nav interface template token = N. Note: when NavUnits in the file header is 0, values are in meters (northings and eastings).  When NavUnits is 3, values are in Lat/Long.  Also see the Layback value, below.
		Head.HSonarStatus(n)=fread(fId,1,'uint16');				%Tow Cable information. System status value, sonar dependant (displayed in Status window).
		Head.HRangeToFish(n)=fread(fId,1,'uint16');				%Slant range to sensor in decimeters; nav interface template token = ? (question mark).  Stored only – not used in any computation.
		Head.HBearingToFish(n)=fread(fId,1,'uint16');			%Bearing to towfish from ship, stored in degrees multiplied by 100; nav interface template token = > (greater-than sign).  Stored only – not used in any computation in Isis.
		%fseek(fId, 2, 0);	%Head.HCableOut(n)=fread(fId,1,'uint16');		%Tow Cable information. Amount of cable payed out in meters; nav interface template token = o.
		%fseek(fId, 4, 0);	%Head.HLayback(n)=fread(fId,1,'float32');		%Tow Cable information. Distance over ground from ship to fish.; nav interface template token = l. Isis Note: When this value is non-zero, Isis assumes that SensorYcoordinate and SensorXcoordinate need to be adjusted with the Layback.  The sensor position is then computed using the current sensor heading and this layback value.  The result is displayed when a position is computed in Isis.
		%fseek(fId, 4, 0);	%Head.HCableTension(n)=fread(fId,1,'float32');	%Tow Cable information Cable tension from serial port. Stored only; nav interface template token = P
		fseek(fId, 10, 0);
		Head.HSensorDepth(n)=fread(fId,1,'float32');			%Sensor Attitude information. Distance (m) from sea surface to sensor. The deeper the sensor goes, the bigger (positive) this value becomes. nav interface template token = 0 (zero)
		Head.HSensorPrimaryAltitude(n)=fread(fId,1,'float32');	%Sensor Attitude information. Distance from towfish to the sea floor; nav interface template token = 7. Isis Note: This is the primary altitude as tracked by the Isis bottom tracker or entered manually by the user. Although not recommended, the user can override the Isis bottom tracker by sending the primary altitude over the serial port.  The user should turn the Isis bottom tracker Off when this is done.
		Head.HSensorAuxAltitude(n)=fread(fId,1,'float32');		%Sensor Attitude information. Auxiliary altitude; nav interface template token = a. Isis Note: This is an auxiliary altitude as transmitted by an altimeter and received over a serial port. The user can switch between the Primary and  Aux altitudes via the "options" button in the Isis bottom track window.
		Head.HSensorPitch(n)=fread(fId,1,'float32');			%Sensor Attitude information. Pitch in degrees (positive=nose up); nav interface template token = 8.
		Head.HSensorRoll(n)=fread(fId,1,'float32');				%Sensor Attitude information. Roll in degrees (positive=roll to starboard); nav interface template token = 9.
		Head.HSensorHeading(n)=fread(fId,1,'float32');			%Sensor Attitude information. Sensor heading in degrees; nav interface template token = h.
		Head.HHeave(n)=fread(fId,1,'float32');					%Attitude information. Sensors heave at start of ping. Positive value means sensor moved up. Note: These Pitch, Roll, Heading, Heave and Yaw values are those received closest in time to this sonar or bathymetry update.  If a TSS or MRU is being used with a multibeam/bathymetry sensor, the user should  use the higher-resolution attitude data found in the XTFATTITUDEDATA structures.
		Head.HYaw(n)=fread(fId,1,'float32');					%Attitude information. Sensor yaw.  Positive means turn to right. Note: These Pitch, Roll, Heading, Heave and Yaw values are those received closest in time to this sonar or bathymetry update.  If a TSS or MRU is being used with a multibeam/bathymetry sensor, the user should use the higher-resolution attitude data found in the XTFATTITUDEDATA structures.  Since the heading information is updated in high resolution, it is not necessary to log or use Yaw in any processing.  Isis does not use Yaw
		%fseek(fId, 4, 0);	%Head.HAttitudeTimeTag(n)=fread(fId,1,'uint32');		%Attitude information. In milliseconds - used to coordinate with millisecond time value in Attitude packets.  (M)andatory when logging XTFATTITUDE packets.
		%fseek(fId, 4, 0);	%Head.HDOT(n)=fread(fId,1,'float32');					%Misc. Distance Off Track
		%fseek(fId, 4, 0);	%Head.HNavFixMilliseconds(n)=fread(fId,1,'uint32');		%Misc. millisecond clock value when nav received
		%fseek(fId, 1, 0);	%Head.HComputerClockHour(n)=fread(fId,1,'uint8');		%Isis Note: The Isis computer clock time when this ping was received. May be different from ping time at start of this record if the sonar time-stamped the data and the two systems aren't synched. This time should be ignored in most cases
		%fseek(fId, 1, 0);	%Head.HComputerClockMinute(n)=fread(fId,1,'uint8');		%Isis Note:  see above Isis Note
		%fseek(fId, 1, 0);	%Head.HComputerClockSecond(n)=fread(fId,1,'uint8');		%Isis Note:  see above Isis Note
		%fseek(fId, 1, 0);	%Head.HComputerClockHsec(n)=fread(fId,1,'uint8');		%Isis Note:  see above Isis Note
		%fseek(fId, 2, 0);	%Head.HFishPositionDeltaX(n)=fread(fId,1,'int16');		%Additional Tow Cable and Fish information from Trackpoint. Stored as meters multiplied by 3.0, supporting +/- 10000.0m (usually from trackpoint); nav interface template token = {DX}.
		%fseek(fId, 2, 0);	%Head.HFishPositionDeltaY(n)=fread(fId,1,'int16');		%Additional Tow Cable and Fish information from Trackpoint. X, Y offsets can be used instead of logged layback.; nav interface template token = {DY}.
		%fseek(fId, 1, 0);	%Head.HFishPositionErrorCode(n)=fread(fId,1,'*char');	%Additional Tow Cable and Fish information from Trackpoint. Error code for FishPosition delta x,y. (typically reported by Trackpoint).
		%fseek(fId, 4, 0);	%Head.HOptionalOffsey(n)=fread(fId,1,'uint32');			%OptionalOffsey (Triton 7125 only)
		%fseek(fId, 1, 0);	%Head.HCableOutHundredths(n)=fread(fId,1,'uint8');		%Hundredths of a meter of cable out, to be added  to the CableOut field.
		%fseek(fId, 6, 0);	%Head.HReservedSpace2(:,n)=fread(fId,6,'uint8');		%Unused. Set to 0.
		fseek(fId, 32, 0);

		for (nn = 1:XtfHead.RNumChansToFollow(nHead(n)))		% Begin ChanInfo Read
			Head.CChannelNumber(nn,n)=fread(fId,1,'uint16');	%Typically 0=port (low frequency) 1=stbd (low frequency) 2=port (high frequency) 3=stbd (high frequency)
			fseek(fId, 2, 0);	%Head.CDownsampleMethod(nn,n)=fread(fId,1,'uint16');	%2 = MAX; 4 = RMS
			Head.CSlantRange(nn,n)=fread(fId,1,'float32');		%Slant range of the data in meters
			Head.CGroundRange(nn,n)=fread(fId,1,'float32');		%Ground range of the data; in meters (SlantRange^2 - Altitude^2)
			%fseek(fId, 4, 0);	%Head.CTimeDelay(nn,n)=fread(fId,1,'float32');		%Amount of time, in seconds, to the start of recorded data. (almost always 0.0).
			%fseek(fId, 4, 0);	%Head.CTimeDuration(nn,n)=fread(fId,1,'float32');	%Amount of time, in seconds, recorded (typically SlantRange/750)
			%fseek(fId, 4, 0);	%Head.CSecondsPerPing(nn,n)=fread(fId,1,'float32');	%Amount of time, in seconds, from ping to ping. (SlantRange/750)
			%fseek(fId, 2, 0);	%Head.CProcessingFlags(nn,n)=fread(fId,1,'uint16');	%4 = TVG; 8 = BAC&GAC; 16 = filter, etc. (almost always zero)
			%fseek(fId, 2, 0);	%Head.CFrequency(nn,n)=fread(fId,1,'uint16');		%Ccenter transmit frequency for this channel
			%fseek(fId, 2, 0);	%Head.CInitialGainCode(nn,n)=fread(fId,1,'uint16');	%Settings as transmitted by sonar
			%fseek(fId, 2, 0);	%Head.CGainCode(nn,n)=fread(fId,1,'uint16');			%Settings as transmitted by sonar
			%fseek(fId, 2, 0);	%Head.CBandWidth(nn,n)=fread(fId,1,'uint16');		%Settings as transmitted by sonar
			%fseek(fId, 4, 0);	%Head.CContactNumber(nn,n)=fread(fId,1,'uint32');	%Contact information . Upated when contacts are saved in Target utility.
			%fseek(fId, 2, 0);	%Head.CContactClassification(nn,n)=fread(fId,1,'uint16'); %Contact information . Updated when contacts are saved in Target utility.
			%fseek(fId, 1, 0);	%Head.CContactSubNumber(nn,n)=fread(fId,1,'uint8');	%Contact information . Udated when contacts are saved in Target utility
			%fseek(fId, 1, 0);	%Head.CContactType(nn,n)=fread(fId,1,'uint8');		%Contact information . Updated when contacts are saved in Target utility
			fseek(fId, 30, 0);
			Head.CNumSamples(nn,n)=fread(fId,1,'uint32');		%Number of samples that will follow this structure. The number of bytes will be this value multiplied by the number of bytes per sample. BytesPerSample found in CHANINFO structure (given in the file header).
			%fseek(fId, 2, 0);	%Head.CMillivoltScale(nn,n)=fread(fId,1,'uint16');	%Maximum voltage, in mv, represented by a full-scale value in the data.If zero, then the value stored in the VoltScale should be used instead. VoltScale can be found in the XTF file header, ChanInfo structure. Note that VoltScale is specified in volts, while MillivoltScale is stored in millivolts. This provides for a range of –65,536 volts to 65,535 volts.
			%fseek(fId, 4, 0);	%Head.CContactTimeOffTrack(nn,n)=fread(fId,1,'float32'); %Time off track to this contact (stored in milliseconds)
			%fseek(fId, 1, 0);	%Head.CContactCloseNumber(nn,n)=fread(fId,1,'uint8');
			%fseek(fId, 1, 0);	%Head.CReserved2(nn,n)=fread(fId,1,'uint8');		%Unused. Set to 0.
			%fseek(fId, 4, 0);	%Head.CFixedVSOP(nn,n)=fread(fId,1,'float32');		%This is the fixed, along-track size of each ping, stored in centimeters. On multibeam systems with zero beam spread, this value needs to be filled in to prevent Isis from calculating along-track ground coverage based on beam spread and speed over ground.
			%fseek(fId, 2, 0);	%Head.CWeight(nn,n)=fread(fId,1,'int16');			%Weighting factor passed by some sonars, this value is mandatory for Edgetech digital sonars types 24, 35, 38, 48 and Kongsberg SA type 48
			%fseek(fId, 4, 0);	%Head.CReservedSpace(nn,n,:)=fread(fId,4,'uint8');	%Unused. Set to 0.
			fseek(fId, 18, 0);
			% Begin Data pass
			CDataSeek(nn,n) = ftell(fId);
			fseek(fId, Head.CNumSamples(nn,n) .* XtfHead.CBytesPerSample(nn), 'cof');
		end
		df = ftell(fId);
	end

	%Data = nan(max(Head.CNumSamples(:)), LenHead, ChF);			% Pre-allocate
	L = (find(XtfHead.RSubChannelNumber(1) == XtfHead.CSubChannelNumber));
	DBit = XtfHead.Descript.BytesPerSample.C{XtfHead.Descript.UniPolar.Code == XtfHead.CUniPolar(L(1)), ...
	                                         XtfHead.Descript.BytesPerSample.Code == XtfHead.CBytesPerSample(L(1))};
	Data = alloc_mex(max(Head.CNumSamples(:)), LenHead, ChF, DBit);			% Pre-allocate

	df=0;	fseek(fId, 0, 'bof');
	for (n = 1:LenHead)			% Begin Data Read
		for (nn = 1:XtfHead.RNumChansToFollow(nHead(n)))
			fseek(fId,CDataSeek(nn,n)-df, 'cof');
			%L = (find(XtfHead.RSubChannelNumber(n) == XtfHead.CSubChannelNumber));
			%DBit = XtfHead.Descript.BytesPerSample.C{XtfHead.Descript.UniPolar.Code == XtfHead.CUniPolar(L(1)), XtfHead.Descript.BytesPerSample.Code == XtfHead.CBytesPerSample(L(1))};
			Data(1:Head.CNumSamples(nn,n),n,nn) = fread(fId, Head.CNumSamples(nn,n), DBit);
			df = ftell(fId);
		end
	end
	fclose(fId);
