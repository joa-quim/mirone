function XtfHead=gFXtfReadHeader(fName,flStat)
%Read XtfHead structure (include XTFFILEHEADER, CHANINFO and Short Header) from *.xtf file.
%
% XtfHead - XtfHead structure;
% fName - the target file name;
% flStat - flag for statistics display (1 or 0).
% Short Header include the first six fields, general for most Headers: XtfHead.RFace, XtfHead.RHeaderType, XtfHead.RSubChannelNumber, XtfHead.RNumChansToFollow, XtfHead.RUnused, XtfHead.RNumBytesThisRecord.
% XtfHead include the next addition fields: XtfHead.Descript, XtfHead.fName, XtfHead.RSeek, XtfHead.ROnFlag.
% Example: XtfHead=gFXtfReadHeader('c:\temp\1.xtf',1);
%
% mail@ge0mlib.ru 01/08/2016
% This function comes from geomlib (FEX: 58387, BSD Licence)

% $Id: gFXtfReadHeader.m 10387 2018-04-27 16:09:09Z j $

	Descript.SonarType.Code=(0:58)';
	Descript.SonarType.Text={'0=NONE , defaul','1=JAMSTEC, Jamstec chirp 2-channel subbottom','2=ANALOG_C31, PC31 8-channel','3=SIS1000, Chirp SIS-1000 sonar','4=ANALOG_32CHAN, Spectrum with 32-channel DSPlink card','5=KLEIN2000, Klein system 2000 with digital interface','6=RWS, Standard PC31 analog with special nav code','7=DF1000, EG&G DF1000 digital interface','8=SEABAT, Reson SEABAT 900x analog/serial','9=KLEIN595, 4-chan Klein 595, same  as ANALOG_C31',...
    '10=EGG260, 2-channel EGG260, same as ANALOG_C31','11=SONATECH_DDS, Sonatech Diver Detection System on Spectrum DSP32C','12=ECHOSCAN, Odom EchoScanII multibeam (with simultaneous analog sidescan)','13=ELAC, Elac multibeam system','14=KLEIN5000, Klein system 5000 with digital interface','15=Reson Seabat 8101','16=Imagenex model 858','17=USN SILOS with 3-channel analog','18=Sonatech Super-high res sidescan sonar','19=Delph AU32 Analog input (2 channel)',...
    '20=Generic sonar using the memory-mapped file interface','21=Simrad SM2000 Multibeam Echo Sounder','22=Standard multimedia audio','23=Edgetech (EG&G) ACI card for 260 sonar through PC31 card','24=Edgetech Black Box','25=Fugro deeptow','26=C&C''s Edgetech Chirp conversion program','27=DTI SAS Synthetic Aperture processor (memmap file)','28=Fugro''s Osiris AUV Sidescan data','29=Fugro''s Osiris AUV Multibeam data',...
    '30=Geoacoustics SLS','31=Simrad EM2000/EM3000','32=Klein system 3000','33=SHRSSS Chirp system','34=Benthos C3D SARA/CAATI','35=Edgetech MP-X','36=CMAX','37=Benthos sis1624','38=Edgetech 4200','39=Benthos SIS1500',...
    '40=Benthos SIS1502','41=Benthos SIS3000','42=Benthos SIS7000','43=DF1000 DCU','44=NONE_SIDESCAN','45=NONE_MULTIBEAM','46=Reson 7125','47=CODA Echoscope','48=Kongsberg SAS','49=QINSy',...
    '50=GeoAcoustics  DSSS','51=CMAX_USB','52=SwathPlus Bathy','53=R2Sonic QINSy','54=R2Sonic Triton','55=Converted SwathPlus Bathy','56=Edgetech 4600','57=Klein 3500','58=Klein 5900'}';

	Descript.NavUnits.Code=[0,3]';
	Descript.NavUnits.Text={'0=Meters (i.e., UTM)','3=Lat/Long'}';
	Descript.TypeOfChannel.Code=(0:3)';
	Descript.TypeOfChannel.Text={'0=SUBBOTTOM','1=PORT','2=STBD','3=BATHYMETRY'}';
	Descript.CorrectionFlags.Code=[1,2]';
	Descript.CorrectionFlags.Text={'1=sonar imagery stored as slant-range','2=sonar imagery stored as ground range (corrected)'}';
	Descript.UniPolar.Code=[0,1]';
	Descript.UniPolar.Text={'0=data is polar','1=data is unipolar'}';
	Descript.BytesPerSample.Code=[1,2,4]';
	Descript.BytesPerSample.Text={'1= 8-bit data','2= 16-bit data','4= 32-bit data'}';
	Descript.BytesPerSample.C={'int8','int16','int32';'uint8','uint16','uint32'};
	Descript.HeaderType.Code=[0:9,10:19,20:28,50,60,61,62,65:69,70:73,100,102:108,199,200]';

	Descript.HeaderType.Text={'0=XTF_HEADER_SONAR (Sidescan data)','1=XTF_HEADER_NOTES','2=XTF_HEADER_BATHY (bathymetry data)','3=XTF_HEADER_ATTITUDE (attitude packet)','4=XTF_HEADER_FORWARD Forward look data (Sonatech)','5=XTF_HEADER_ELAC Elac raw data packet','6=XTF_HEADER_RAW_SERIAL Raw ASCII serial port data','7=XTF_HEADER_EMBED_HEAD Embedded header record - num samples probably changed','8=XTF_HEADER_HIDDEN_SONAR Redundant (overlapping) ping from Klein 5000','9=XTF_HEADER_SEAVIEW_PROCESSED_BATHY Bathymetry (angles) for Seaview',...
    '10=XTF_HEADER_SEAVIEW_DEPTHS Bathymetry from Seaview data (depths)','11=XTF_HEADER_RSVD_HIGHSPEED_SENSOR Used by Klein. 0=roll, 1=yaw.','12=XTF_HEADER_ECHOSTRENGTH Elac EchoStrength (10 values)','13=XTF_HEADER_GEOREC Used to store mosaic parameters','14=XTF_HEADER_KLEIN_RAW_BATHY Bathymetry data from the Klein 5000','15=XTF_HEADER_HIGHSPEED_SENSOR2 High speed sensor from Klein 5000','16=XTF_HEADER_ELAC_XSE Elac dual-head','17=XTF_HEADER_BATHY_XYZA','18=XTF_HEADER_K5000_BATHY_IQ Raw IQ data from Klein 5000 server','19=XTF_HEADER_BATHY_SNIPPET',...
    '20=XTF_HEADER_GPS GPS Position','21=XTF_HEADER_STAT GPS statistics','22=XTF_HEADER_SINGLEBEAM','23=XTF_HEADER_GYRO Heading/Speed Sensor','24=XTF_HEADER_TRACKPOINT','25=XTF_HEADER_MULTIBEAM','26=XTF_HEADER_Q_SINGLEBEAM','27=XTF_HEADER_Q_MULTITX','28=XTF_HEADER_Q_MULTIBEAM',...
    '50=XTF_HEADER_TIME',...
    '60=XTF_HEADER_BENTHOS_CAATI_SARA. Custom Benthos data','61=XTF_HEADER_7125 7125 Bathy Data','62=XTF_HEADER_7125_SNIPPET 7125 Bathy Data Snippets','65=XTF_HEADER_QINSY_R2SONIC_BATHY QINSy R2Sonic bathymetry data','66=XTF_HEADER_QINSY_R2SONIC_FTS QINSy R2Sonics Foot Print Time Series (snippets)','67=','68=XTF_HEADER_R2SONIC_BATHY Triton R2Sonic bathymetry data','69=XTF_HEADER_R2SONIC_FTS Triton R2Sonic Footprint Time Series',...
    '70=XTF_HEADER_CODA_ECHOSCOPE_DATA Custom CODA Echoscope Data','71=XTF_HEADER_CODA_ECHOSCOPE_CONFIG Custom CODA Echoscope Data','72=XTF_HEADER_CODA_ECHOSCOPE_IMAGE Custom CODA Echoscope Data','73=XTF_HEADER_EDGETECH_4600',...
    '100=XTF_HEADER_POSITION Raw position packet - Reserved for use by Reson, Inc. RESON ONLY','102=XTF_HEADER_BATHY_PROC','103=XTF_HEADER_ATTITUDE_PROC','104=XTF_HEADER_SINGLEBEAM_PROC','105=XTF_HEADER_AUX_PROC Aux Channel + AuxAltitude + Magnetometer','106=XTF_HEADER_KLEIN3000_DATA_PAGE','107=XTF_HEADER_POS_RAW_NAVIGATION','108=XTF_HEADER_KLEINV4_DATA_PAGE',...
    '199=XTF_RAW_CUSTOM_HEADER custom vendor data follows (defined in Xtf.h)',...
    '200=XTF_HEADER_USERDEFINED This packet type is reserved for specific applications. (defined in Xtf.h)'}';

	[fId, mes] = fopen(fName,'r');
	if ~isempty(mes),	error(['gFXtfReadHeader: ' mes]),	end
	finfo = dir(fName);	fSize =	finfo.bytes;

	%===Begin XTFFILEHEADER Structure Read
	XtfHead.Descript=Descript;
	XtfHead.fName=fName;
	XtfHead.HFileFormat=fread(fId,1,'uint8');
	if (XtfHead.HFileFormat~=123),	error('Error gFXtfReadHead: FileFormat~=7B');	end %Set to 123 (0x7B)
	XtfHead.HSystemType=fread(fId,1,'uint8');		%Set to 1
	XtfHead.HRecordingProgramName=fread(fId,8,'*char')'; %Example: "Isis
	XtfHead.HRecordingProgramVersion=fread(fId,8,'*char')'; %Example: "556" for version 5.56
	XtfHead.HSonarName=fread(fId,16,'*char')';		%Name of server used to access sonar.  Example: "C31_SERV.EXE"
	XtfHead.HSonarType=fread(fId,1,'uint16');
	XtfHead.HNoteString=fread(fId,64,'*char')';		%Notes as entered in the Sonar Setup dialog box
	XtfHead.HThisFileName=fread(fId,64,'*char')';	%Name of this file. Example:"LINE12-B.XTF"
	XtfHead.HNavUnits=fread(fId,1,'uint16');		%0=Meters (i.e., UTM) or 3=Lat/Long
	XtfHead.HNumberOfSonarChannels=fread(fId,1,'uint16'); %if > 6, header grows to 2K in size
	XtfHead.HNumberOfBathymetryChannels=fread(fId,1,'uint16');
	XtfHead.HNumberOfSnippetChannels=fread(fId,1,'uint8');
	XtfHead.HNumberOfForwardLookArrays=fread(fId,1,'uint8');
	XtfHead.HNumberOfEchoStrengthChannels=fread(fId,1,'uint16');
	XtfHead.HNumberOfInterferometryChannels=fread(fId,1,'uint8');
	XtfHead.HReserved1=fread(fId,1,'uint8');		%Reserved. Set to 0.
	XtfHead.HReserved2=fread(fId,1,'uint16');		%Reserved. Set to 0.
	XtfHead.HReferencePointHeigh=fread(fId,1,'float32'); %Height of reference point above water line (m)
	XtfHead.HProjectionType=fread(fId,12,'*char')';	%Navigation System Parameters. Not currently used. Set to 0.
	XtfHead.HSpheriodType=fread(fId,10,'*char')';	%Navigation System Parameters. Not currently used. Set to 0.
	XtfHead.HNavigationLatency=fread(fId,1,'int32');%Navigation System Parameters. Latency of nav system in milliseconds. (Usually GPS). ISIS Note: This value is entered on the Serial port setup dialog box.  When computing a position, Isis will take the time of the navigation and subtract this value.
	XtfHead.HOriginY=fread(fId,1,'float32');		%Navigation System Parameters. Not currently used. Set to 0
	XtfHead.HOriginX=fread(fId,1,'float32');		%Navigation System Parameters. Not currently used. Set to 0
	XtfHead.HNavOffsetY=fread(fId,1,'float32');		%Navigation System Parameters. Orientation of positive Y is forward. ISIS Note: This offset is entered in the Multibeam setup dialog box
	XtfHead.HNavOffsetX=fread(fId,1,'float32');		%Navigation System Parameters. Orientation of positive X is to starboard. ISIS Note: This offset is entered in the Multibeam setup dialog box
	XtfHead.HNavOffsetZ=fread(fId,1,'float32');		%Navigation System Parameters. Orientation of positive Z is down.  Just like depth. ISIS Note: This offset is entered in the Multibeam setup dialog box
	XtfHead.HNavOffsetYaw=fread(fId,1,'float32');	%Navigation System Parameters. Orientation of positive yaw is turn to right. ISIS Note: This offset is entered in the Multibeam setup dialog box
	XtfHead.HMRUOffsetY=fread(fId,1,'float32');		%Navigation System Parameters. Orientation of positive Y is forward. ISIS Note: This offset is entered in the Multibeam setup dialog box
	XtfHead.HMRUOffsetX=fread(fId,1,'float32');		%Navigation System Parameters. Orientation of positive X is to starboard. ISIS Note: This offset is entered in the Multibeam setup dialog box
	XtfHead.HMRUOffsetZ=fread(fId,1,'float32');		%Navigation System Parameters. Orientation of positive Z is down.  Just like depth. ISIS Note: This offset is entered in the Multibeam setup dialog box
	XtfHead.HMRUOffsetYaw=fread(fId,1,'float32');	%Navigation System Parameters. Orientation of positive yaw is turn to right. ISIS Note: This offset is entered in the Multibeam setup dialog box
	XtfHead.HMRUOffsetPitch=fread(fId,1,'float32');	%Navigation System Parameters. Orientation of positive pitch is nose up. ISIS Note: This offset is entered in the Multibeam setup dialog box. ISIS Note: This offset is entered in the Multibeam setup dialog box
	XtfHead.HMRUOffsetRoll=fread(fId,1,'float32');	%Navigation System Parameters. Orientation of positive roll is lean to starboard. ISIS Note: This offset is entered in the Multibeam setup dialog box
	nChanel = XtfHead.HNumberOfSonarChannels + XtfHead.HNumberOfBathymetryChannels + XtfHead.HNumberOfSnippetChannels + ...
	          XtfHead.HNumberOfForwardLookArrays + XtfHead.HNumberOfEchoStrengthChannels + XtfHead.HNumberOfInterferometryChannels;

	% Begin CHANINFO Structure Allocate
	v = zeros(1, nChanel);
	XtfHead.CTypeOfChannel = v;
	XtfHead.CSubChannelNumber = v;
	XtfHead.CCorrectionFlags = v;
	XtfHead.CUniPolar = v;
	XtfHead.CBytesPerSample = v;
	XtfHead.CReserved = v;
	XtfHead.CChannelName=char(zeros(16,nChanel));
	XtfHead.CVoltScale = v;
	XtfHead.CFrequency = v;
	XtfHead.CHorizBeamAngle = v;
	XtfHead.CTiltAngle = v;
	XtfHead.CBeamWidth = v;
	XtfHead.COffsetX = v;
	XtfHead.COffsetY = v;
	XtfHead.COffsetZ = v;
	XtfHead.COffsetYaw = v;
	XtfHead.COffsetPitch = v;
	XtfHead.COffsetRoll = v;
	XtfHead.CBeamsPerArray = v;
	XtfHead.CReservedArea2=char(zeros(54,nChanel));

	for (m = 1:nChanel)		% Begin CHANINFO Structure Read
		XtfHead.CTypeOfChannel(m)=fread(fId,1,'uint8'); %SUBBOTTOM=0, PORT=1, STBD=2, BATHYMETRY=3
		XtfHead.CSubChannelNumber(m)=fread(fId,1,'uint8'); %Index for which CHANINFO structure this is
		XtfHead.CCorrectionFlags(m)=fread(fId,1,'uint16'); %1=sonar imagery stored as slant-range, 2=sonar imagery stored as ground range (corrected)
		XtfHead.CUniPolar(m)=fread(fId,1,'uint16'); %0=data is polar, 1=data is unipolar
		XtfHead.CBytesPerSample(m)=fread(fId,1,'uint16'); %1 (8-bit data) or 2 (16-bit data) or 4 (32-bit)
		XtfHead.CReserved(m)=fread(fId,1,'uint32'); %Isis Note: Previously this was SamplesPerChannel.  Isis now supports the recording of every sample per ping, which means that number of samples per channel can vary from ping to ping if the range scale changes.  Because of this, the NumSamples value in the XTFPINGCHANHEADER structure (defined in Section 3.18) holds the number of samples to read for a given channel. For standard analog systems, this Reserved value is still filled in with 1024, 2048 or whatever the initial value is for SamplesPerChannel.
		XtfHead.CChannelName(:,m)=fread(fId,16,'*char')'; %Text describing channel.  i.e., "Port 500"
		XtfHead.CVoltScale(m)=fread(fId,1,'float32'); %This states how many volts are represented by a maximum sample value in the range  [-5.0 to +4.9998] volts. Default is 5.0.
		XtfHead.CFrequency(m)=fread(fId,1,'float32'); %Center transmit frequency
		XtfHead.CHorizBeamAngle(m)=fread(fId,1,'float32'); %Typically 1 degree or so
		XtfHead.CTiltAngle(m)=fread(fId,1,'float32'); %Typically 30 degrees
		XtfHead.CBeamWidth(m)=fread(fId,1,'float32'); %3dB beam width, Typically 50 degrees
		XtfHead.COffsetX(m)=fread(fId,1,'float32'); %Orientation of positive X is to starboard. Note: This offset is entered in the Multibeam setup dialog box
		XtfHead.COffsetY(m)=fread(fId,1,'float32'); %Orientation of positive Y is forward. Note: This offset is entered in the Multibeam setup dialog box
		XtfHead.COffsetZ(m)=fread(fId,1,'float32'); %Orientation of positive Z is down.  Just like depth. Note: This offset is entered in the Multibeam setup dialog box
		XtfHead.COffsetYaw(m)=fread(fId,1,'float32'); %Orientation of positive yaw is turn to right. If the multibeam sensor is reverse mounted (facing backwards), then OffsetYaw will be around 180 degrees. Note: This offset is entered in the Multibeam setup dialog box
		XtfHead.COffsetPitch(m)=fread(fId,1,'float32'); %Orientation of positive pitch is nose up. Note: This offset is entered in the Multibeam setup dialog box
		XtfHead.COffsetRoll(m)=fread(fId,1,'float32'); %Orientation of positive roll is lean to starboard. Note: This offset is entered in the Multibeam setup dialog box
		XtfHead.CBeamsPerArray(m)=fread(fId,1,'uint16'); %For forward look only (i.e., Sonatech DDS)
		XtfHead.CReservedArea2(:,m)=fread(fId,54,'*char')'; %Unused Set value to 0
	end

	fseek(fId,1024-mod(ftell(fId),1024),'cof');
	ss=ftell(fId);nRec=0;
	while (fSize > ftell(fId))		% Begin Calc Num of ShortHeader structure
		nRec=nRec+1;
		face=fread(fId,1,'uint16');
		if (face ~= 64206),		error('Error gFXtfReadHead: MagicNumber~=FACE');	end
		fseek(fId,8,'cof');
		HNumBytesThisRecord = fread(fId,1,'uint32'); %Total byte count for this ping including this ping header. Isis Note: Isis records data packets in multiples of 64 bytes. If the size of the data packet is not an exact multiple of 64 bytes, zeros are padded at the end packet and this value will be promoted to the next 64-byte granularity.  In all cases, this value will be the EXACT size of this packet.
		fseek(fId,HNumBytesThisRecord-14,'cof');
	end
	
	%===Begin ShortHeader structure Allocate
	v = zeros(1, nRec);
	XtfHead.RFace = v;
	XtfHead.RHeaderType = v;
	XtfHead.RSubChannelNumber = v;
	XtfHead.RNumChansToFollow = v;
	XtfHead.RUnused = v;
	XtfHead.RNumBytesThisRecord = v;
	XtfHead.RSeek = v;
	XtfHead.ROnFlag = v;

	fseek(fId,ss,'bof');
	for (m = 1:nRec)		% Begin ShortHeader structure Read
		XtfHead.RFace=fread(fId,1,'uint16');if XtfHead.RFace~=64206, error('Error gFXtfReadHead: MagicNumber~=FACE');	end
		XtfHead.RHeaderType(m)=fread(fId,1,'uint8');
		XtfHead.RSubChannelNumber(m)=fread(fId,1,'uint8'); %If HeaderType is bathymetry, this indicates which head; if HeaderType is forward-looking sonar, and then this indicates which array. Also, Klein 5000 beam numbers are logged here. 
		XtfHead.RNumChansToFollow(m)=fread(fId,1,'uint16'); %If HeaderType is sonar, number of channels to follow.
		XtfHead.RUnused(m)=fread(fId,1,'uint32'); %Unused. Set to 0
		XtfHead.RNumBytesThisRecord(m)=fread(fId,1,'uint32'); %Total byte count for this ping including this ping header. Isis Note: Isis records data packets in multiples of 64 bytes. If the size of the data packet is not an exact multiple of 64 bytes, zeros are padded at the end packet and this value will be promoted to the next 64-byte granularity.  In all cases, this value will be the EXACT size of this packet.
		XtfHead.ROnFlag(m)=1;
		XtfHead.RSeek(m)=ftell(fId);
		fseek(fId,XtfHead.RNumBytesThisRecord(m)-14,'cof');
	end
	fclose(fId);

	%===Begin Statistics display
	if flStat
		a1=sparse(XtfHead.RHeaderType+1,ones(1,nRec),ones(1,nRec));a1Mess=find(a1)-1;a1Num=nonzeros(a1);
		for n=1:size(a1Mess,1)
			if any(XtfHead.Descript.HeaderType.Code==a1Mess(n))
				s = XtfHead.Descript.HeaderType.Text{XtfHead.Descript.HeaderType.Code == a1Mess(n)};
			else
				s = [num2str(a1Mess(n)) '=Not Defined'];
			end
			fprintf('Rec: %s, Num: %d ',s,a1Num(n));
			L=find(XtfHead.RHeaderType==a1Mess(n));L2=size(L,2);
			a2=sparse(XtfHead.RSubChannelNumber(L)+1,ones(1,L2),ones(1,L2));a2Mess=find(a2)-1;a2Num=nonzeros(a2);
			for (nn = 1:size(a2Mess,1))
				fprintf('[ SubCh: %d, Num: %d;',a2Mess(nn),a2Num(nn));
				L=find((XtfHead.RHeaderType==a1Mess(n))&(XtfHead.RSubChannelNumber==a2Mess(nn)));L2=size(L,2);
				a3=sparse(XtfHead.RNumChansToFollow(L)+1,ones(1,L2),ones(1,L2));a3Mess=find(a3)-1;a3Num=nonzeros(a3);
				fprintf(' ChFollow:');	fprintf(' %d',a3Mess);
				fprintf(', Num:');		fprintf(' %d',a3Num);	fprintf(' ]');
			end
			fprintf('\n');
		end
	end
