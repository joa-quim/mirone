%#mex
function  varargout = hdf_funs(fun,varargin)
% ...
	[varargout{1:nargout}] = feval(fun, varargin{:});

% ---------------------------------------------------------------------
function fileinfo = hdfinfo(varargin)
%HDFINFO Information about HDF 4 or HDF-EOS 2 file
%
%	A SELF-CONTAINED VERSION OF HDFINFO THAT RELIES ONLY ON HDF.DLL -- JL 2009
%
%   FILEINFO = HDFINFO(FILENAME) returns a structure whose fields contain
%   information about the contents an HDF or HDF-EOS file.  FILENAME is a
%   string that specifies the name of the HDF file. HDF-EOS files are
%   described as HDF files.
%
%   FILEINFO = HDFINFO(FILENAME,MODE) reads the file as an HDF file if MODE
%   is 'hdf', or as an HDF-EOS file if MODE is 'eos'.  If MODE is 'eos',
%   only HDF-EOS data objects are queried.  To retrieve information on the
%   entire contents of a hybrid HDF-EOS file, MODE must be 'hdf' (default).
%   
%   The set of fields in FILEINFO depends on the individual file.  Fields
%   that may be present in the FILEINFO structure are:
%
%   HDF objects:
%   
%   Filename   A string containing the name of the file
%   Vgroup     An array of structures describing the Vgroups
%   SDS        An array of structures describing the Scientific Data Sets
%   Vdata      An array of structures describing the Vdata sets
%   Raster8    An array of structures describing the 8-bit Raster Images
%   Raster24   An array of structures describing the 24-bit Raster Images
%   
%   HDF-EOS objects:
%
%   Point      An array of structures describing HDF-EOS Point data
%   Grid       An array of structures describing HDF-EOS Grid data
%   Swath      An array of structures describing HDF-EOS Swath data
%   
%   The data set structures above share some common fields.  They are (note,
%   not all structures will have all these fields):
%   
%   Filename          A string containing the name of the file
%   Type              A string describing the type of HDF object 
%   Name              A string containing the name of the data set
%   Attributes        An array of structures with fields 'Name' and 'Value'
%                     describing the name and value of the attributes of the data set
%   Rank              A number specifying the number of dimensions of the data set
%   Ref               The reference number of the data set
%   Label             A cell array containing an Annotation label
%   Description       A cell array containing an Annotation description
%
%   Fields specific to each structure are:
%   
%   Vgroup:
%   
%      Class      A string containing the class name of the data set
%      Vgroup     An array of structures describing Vgroups
%      SDS        An array of structures describing Scientific Data sets
%      Vdata      An array of structures describing Vdata sets
%      Raster24   An array of structures describing 24-bit raster images  
%      Raster8    An array of structures describing 8-bit raster images
%      Tag        The tag of this Vgroup
%                 
%   SDS:
%              
%      Dims       An array of structures with fields 'Name', 'DataType',
%                 'Size', 'Scale', and 'Attributes'.  Describing the
%                 dimensions of the data set.  'Scale' is an array of numbers
%                 to place along the dimension and demarcate intervals in the data set.
%      DataType   A string specifying the precision of the data
%      Index      Number indicating the index of the SDS
%   
%   Vdata:
%   
%      DataAttributes    An array of structures with fields 'Name' and 'Value'
%                        describing the name and value of the attributes of the entire data set
%      Class             A string containing the class name of the data set
%      Fields            An array of structures with fields 'Name' and
%                        'Attributes' describing the fields of the Vdata
%      NumRecords        A number specifying the number of records of the data set   
%      IsAttribute       1 if the Vdata is an attribute, 0 otherwise
%      
%   Raster8 and Raster24:
%
%      Name           A string containing the name of the image
%      Width          An integer indicating the width of the image in pixels
%      Height         An integer indicating the height of the image in pixels
%      HasPalette     1 if the image has an associated palette, 0 otherwise (8-bit only)
%      Interlace      A string describing the interlace mode of the image (24-bit only)
%
%   Point:
%
%      Level          A structure with fields 'Name', 'NumRecords',
%                     'FieldNames', 'DataType' and 'Index'.  This structure
%                     describes each level of the Point
%      
%   Grid:
%     
%      UpperLeft      A number specifying the upper left corner location in meters
%      LowerRight     A number specifying the lower right corner location in meters
%      Rows           An integer specifying the number of rows in the Grid
%      Columns        An integer specifying the number of columns in the Grid
%      DataFields     An array of structures with fields 'Name', 'Rank', 'Dims',
%                     'NumberType', 'FillValue', and 'TileDims'. Each structure
%                     describes a data field in the Grid fields in the Grid
%      Projection     A structure with fields 'ProjCode', 'ZoneCode',
%                     'SphereCode', and 'ProjParam' describing the Projection
%                     Code, Zone Code, Sphere Code and projection parameters of
%                     the Grid
%      Origin Code    A number specifying the origin code for the Grid
%      PixRegCode     A number specifying the pixel registration code
%      
%   Swath:
%		       
%      DataFields         An array of structures with fields 'Name', 'Rank', 'Dims',
%                         'NumberType', and 'FillValue'.  Each structure
%                         describes a Data field in the Swath 
%      GeolocationFields  An array of structures with fields 'Name', 'Rank', 'Dims',
%                         'NumberType', and 'FillValue'.  Each structure
%                         describes a Geolocation field in the Swath 
%      MapInfo            A structure with fields 'Map', 'Offset', and
%                         'Increment' describing the relationship between the
%                         data and geolocation fields. 
%      IdxMapInfo         A structure with 'Map' and 'Size' describing the
%                         relationship between the indexed elements of the
%                         geolocation mapping
%   
% 
%   Example:  
%             % Retrieve info about example.hdf
%             fileinfo = hdfinfo('example.hdf');
%             % Retrieve info about Scientific Data Set in example
%             data_set_info = fileinfo.SDS;
  
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2005/11/15 01:10:11 $

[filename, mode] = parseInputs(varargin{:});

%  Open the interfaces. The private functions will expect the interfaces to be
%  open and closed by this function.   This is done for performance; opening
%  and closing these interfaces is expensive.    
if strcmpi(mode,'hdf')  
  fileID = hdf('H','open',filename,'read',0);
  if fileID==-1
    error('MATLAB:hdfinfo:fileOpen', 'Unable to open file.');
  end
  sdID = hdf('SD','start',filename,'read');
  if sdID==-1
    warning('MATLAB:hdfinfo:sdInterfaceStart', ...
            'Unable to start SD interface. Scientific Data Set information is unavailable.');
  end
  vstatus = hdf('V','start',fileID);
  if vstatus==-1
    warning('MATLAB:hdfinfo:vInterfaceStart', ...
            'Unable to start V interface.  Vgroup and Vdata information is unavailable.');
  end
  anID = hdf('AN','start',fileID);
  if anID == -1
    warning('MATLAB:hdfinfo:anInterfaceStart', ...
            'Unable to start AN interface.  Annotation information is unavailable.');
  end
elseif strcmpi(mode,'eos')
  gdID = hdf('GD','open',filename,'read');
  if gdID==-1
    warning('MATLAB:hdfinfo:gdInterfaceStart', ...
            'Unable to start GD interface.  Grid information is unavailable.');
  end
  swID = hdf('SW','open',filename,'read');
  if swID==-1
    warning('MATLAB:hdfinfo:swInterfaceStart', ...
            'Unable to start SW interface.  Swath information is unavailable.');
  end
  ptID = hdf('PT','open',filename,'read');
  if ptID==-1
    warning('MATLAB:hdfinfo:ptInterfaceStart', ...
            'Unable to start PT interface.  Point information is unavailable.');
  end
end

%Get attributes that apply to entire file
fileAttribute = '';
if ~strcmpi(mode,'eos') 
  if sdID ~= -1
    [ndatasets,nGlobalAttrs,status] = hdf('SD','fileinfo',sdID);
    hdfwarn(status)
    if nGlobalAttrs>0
      for i=1:nGlobalAttrs
	[fileAttribute(i).Name,sdsDataType,count,status] = hdf('SD','attrinfo',sdID,i-1);
	hdfwarn(status)
	[fileAttribute(i).Value, status] = hdf('SD','readattr',sdID,i-1);
	hdfwarn(status)
      end
    end
  end
end

% HDF-EOS files are typically "hybrid", meaning they contain HDF-EOS and
% HDF data objects.  It is impossible to distinguish the additional HDF
% objects that may be in an HDF-EOS file from HDF-EOS object.  The only way
% to read these is to look at the entire file as an HDF file.  
if strcmp(mode,'eos')
  if gdID ~= -1
    grid = GridInfo(filename,gdID);
  end
  if swID ~= -1
    swath = SwathInfo(filename,swID);
  end
  if ptID ~= -1
    point = PointInfo(filename,ptID);
  end
else
  %Get Vgroup structures
  if vstatus ~= -1
    [Vgroup,children] = VgroupInfo(filename,fileID,sdID,anID);
  else
    Vgroup = [];
    children = [];
  end
  
  %Get SDS structures
  if sdID ~= -1
    Sds = SdsInfo(filename,sdID,anID,children);
    else
    Sds = [];
  end
  
  %Get Vdata structures
  if vstatus ~= -1
    vdinfo = VdataInfo(filename,fileID,anID);
  else
    vdinfo = [];
  end

  %Get 8-bit image structures
  raster8 = Raster8Info(filename,children,anID);

  %Get 24-bit image structures
  raster24 = Raster24Info(filename,children,anID);

  %Get file annotations
  if anID ~= -1
    [label, desc] = annotationInfo(anID);
  else
    label = {};
    desc = {};
  end
end

%Populate output structure
hinfo.Filename = filename;
if strcmp(mode,'eos')
  if ~isempty(point),	hinfo.Point = point;	end
  if ~isempty(grid),	hinfo.Grid = grid;		end
  if ~isempty(swath),	hinfo.Swath = swath;	end
else
  if ~isempty(fileAttribute),	hinfo.Attributes = fileAttribute;	end
  if ~isempty(Vgroup),	hinfo.Vgroup = Vgroup;	end
  if ~isempty(Sds),		hinfo.SDS = Sds;		end
  if ~isempty(vdinfo),	hinfo.Vdata = vdinfo;	end
  if ~isempty(raster8),	hinfo.Raster8 = raster8;end
  if ~isempty(raster24),hinfo.Raster24 = raster24;  end
  if ~isempty(label),	hinfo.Label = label;	end
  if ~isempty(desc),	hinfo.Description = desc;	end
end

%Close interfaces
if strcmp(mode,'hdf')
  if sdID ~= -1
    status = hdf('SD','end',sdID);
    hdfwarn(status)
  end
  if vstatus ~= -1
    status = hdf('V','end',fileID);
    hdfwarn(status)
  end
  if anID ~= -1
    status = hdf('AN','end',anID);
    hdfwarn(status)
  end
  if fileID ~= -1
    status = hdf('H','close',fileID);
    hdfwarn(status)
  end
elseif strcmp(mode,'eos')  
  if gdID ~= -1
    hdf('GD','close',gdID);
  end
  if swID ~= -1
    hdf('SW','close',swID);
  end
  if ptID ~= -1
    hdf('PT','close',ptID);
  end
end

fileinfo = hinfo;

%===============================================================
function [filename, mode] = parseInputs(varargin)

%Check input arguments
error(nargchk(1,2,nargin));

if ischar(varargin{1})
  filename = varargin{1};
  %Get full path of the file
  fid = fopen(filename);
  if fid~=-1
    filename = fopen(fid);  
    fclose(fid);
  else
    error('MATLAB:hdfinfo:fileNotFound', 'File not found.');
  end
else
  error('MATLAB:hdfinfo:badFilename', 'FILENAME must be a string.');
end

if nargin == 1
  mode = 'hdf';
elseif nargin == 2
  if strcmpi(varargin{2},'eos');
    mode = 'eos';
  elseif strcmpi(varargin{2},'hdf')
    mode = 'hdf';
  else
    error('MATLAB:hdfinfo:badMode', 'Invalid value for MODE.  MODE must be either ''eos'' or ''hdf''.');
  end
end

if ~hdf('H','ishdf',filename)
  error('MATLAB:hdfinfo:invalidFile', '%s is not a valid HDF file.', filename);
end

%================================================================
function [Vgroup, children] = VgroupInfo(filename,fileID,sdID,anID)
Vgroup = [];

%Find top level (lone) Vgroups
[refArray,maxsize] = hdf('V','lone',fileID,0);
[refArray,maxsize] = hdf('V','lone',fileID,maxsize);
hdfwarn(maxsize)

children.Tag = [];
children.Ref = [];
%Get Vgroup structures (including children)
for i=1:maxsize
  [Vgrouptemp, child] = hdfvgroupinfo(filename,refArray(i),fileID, sdID, anID);
  if ~isempty(Vgrouptemp.Filename)
    Vgroup = [Vgroup Vgrouptemp];
  end
  if ~isempty(child.Tag)
    children.Tag = [children.Tag, child.Tag];
    children.Ref = [children.Ref, child.Ref];
  end
end

%================================================================
function vdinfo = VdataInfo(filename,fileID,anID)
vdinfo = [];

[refArray, maxsize] = hdf('VS','lone',fileID,0);
[refArray, maxsize] = hdf('VS','lone',fileID,maxsize);
hdfwarn(maxsize)
for i=1:length(refArray)
  vdtemp = hdfvdatainfo(filename,fileID,anID,refArray(i));
  %Ignore Vdata's that are attributes
  %Ignore Vdata's that are Attr0.0 class, this is consistent with the 
  %NCSA's Java HDF Viewer
  if ~vdtemp.IsAttribute && ~strcmp(vdtemp.Class,'Attr0.0')
    vdinfo = [vdinfo vdtemp];
  end
end

%================================================================
function Sds = SdsInfo(filename,sdID,anID,children)
%Initialize output to empty
Sds = [];

%Get number of data sets in file. SDS index is zero based
[ndatasets,nGlobalAttrs, status] = hdf('SD','fileinfo',sdID);
hdfwarn(status)
for i=1:ndatasets
  sdsID = hdf('SD','select',sdID,i-1);
  if sdID == -1
    warning('MATLAB:hdfinfo:sdRetrieve', 'Unable to retrieve Scientific Data Set information.');
    return
  end
  ref = hdf('SD','idtoref',sdsID);
  hdfwarn(ref)
  status = hdf('SD','endaccess',sdsID);
  hdfwarn(status)
  %Don't add to hinfo structure if it is already part of a Vgroup
  %Assuming that the tag used for SDS is consistent for each file
  if ~(isused(hdf('ML','tagnum','DFTAG_NDG'),ref,children) || ...
       isused(hdf('ML','tagnum','DFTAG_SD'),ref,children) || ...
       isused(hdf('ML','tagnum','DFTAG_SDG'),ref,children) )
    [sds, IsScale] = hdfsdsinfo(filename,sdID,anID,i-1);
    %Ignore dimension scales
    if ~IsScale
      Sds =[Sds sds];    
    end
  end
end

%================================================================
function raster8 = Raster8Info(filename,children,anID)
raster8 = [];

%Get number of 8-bit raster images
nimages = hdf('DFR8','nimages',filename);

for i=1:nimages
  % It wouldn't seem like this call does anything, but it really does.
  [width, height, hasmap, status] = hdf('DFR8','getdims',filename);
  ref = hdf('DFR8','lastref');
  rinfotemp = hdfraster8info(filename,ref,anID);
  if ~(isused(hdf('ML','tagnum','DFTAG_RIG'),ref,children) || ...
       isused(hdf('ML','tagnum','DFTAG_RI'),ref,children) || ...
       isused(hdf('ML','tagnum','DFTAG_CI'),ref,children) || ...
       isused(hdf('ML','tagnum','DFTAG_CI8'),ref,children) || ...
       isused(hdf('ML','tagnum','DFTAG_RI8'),ref,children))
    raster8 = [raster8 rinfotemp];
  end
end

%Restart the DFR8 interface
status = hdf('DFR8','restart');
hdfwarn(status)

%================================================================
function raster24 = Raster24Info(filename,children,anID)
raster24 = [];
%Get number of 24-bit raster images
nimages = hdf('DF24','nimages',filename);

for i=1:nimages
  % It wouldn't seem like this call does anything, but it really does.
  [width, height, hasmap, status] = hdf('DF24','getdims',filename);
  ref = hdf('DF24','lastref');
  rinfotemp = hdfraster24info(filename,ref,anID);
  if ~(isused(hdf('ML','tagnum','DFTAG_RIG'),ref,children) || ...
       isused(hdf('ML','tagnum','DFTAG_RI'),ref,children) || ...
       isused(hdf('ML','tagnum','DFTAG_CI'),ref,children) || ...
       isused(hdf('ML','tagnum','DFTAG_CI8'),ref,children) || ...
       isused(hdf('ML','tagnum','DFTAG_RI8'),ref,children))
       raster24 = [raster24 rinfotemp];
  end
end

%Restart the DF24 interface
status = hdf('DF24','restart');
hdfwarn(status)

%================================================================x
function pinfo = PointInfo(filename,fileID)
pinfo = [];

[numPoints, pointListLong] = hdf('PT','inqpoint',filename);
if numPoints>0
  pointList = parselist(pointListLong);
  for i=1:numPoints
    pinfotemp = hdfpointinfo(filename,fileID,pointList{i});
    pinfo = [pinfo pinfotemp];
  end
end

%================================================================x
function swinfo = SwathInfo(filename,fileID)
swinfo = [];

[numSwaths, swathListLong] = hdf('SW','inqswath',filename);
if numSwaths>0
  swathList = parselist(swathListLong);
  for i=1:numSwaths
    swinfotemp = hdfswathinfo(filename,fileID,swathList{i});
    swinfo = [swinfo swinfotemp];
  end
end

%================================================================x
function gdinfo = GridInfo(filename,fileID)
gdinfo = [];

[numGrids, gridListLong] = hdf('GD','inqgrid',filename);
if numGrids>0
  gridList = parselist(gridListLong);
  for i=1:numGrids
    gdinfotemp = hdfgridinfo(filename,fileID,gridList{i});
    gdinfo = [gdinfo gdinfotemp];
  end
end

%================================================================
function [label,desc] = annotationInfo(anID)
% Retrieves annotations for the file.  
label ={};
desc = {};

[numFileLabel,numFileDesc,numDataLabel,numDataDesc,status] = hdf('AN','fileinfo',anID);
hdfwarn(status)
if status==0
  if numFileLabel>0
    for i=1:numFileLabel
      FileLabelID = hdf('AN','select',anID,i-1,'file_label');
      hdfwarn(FileLabelID)
      if FileLabelID~=-1
	length = hdf('AN','annlen',FileLabelID);
	hdfwarn(length)
	[label{i},status] = hdf('AN','readann',FileLabelID,length);
	hdfwarn(status)
	status = hdf('AN','endaccess',FileLabelID);
	hdfwarn(status) 
      end
    end
  end
  if numFileDesc>0
    for i=1:numFileDesc
      FileDescID = hdf('AN','select',anID,i-1,'file_desc');
      hdfwarn(FileDescID)
      if FileDescID~=-1
	len = hdf('AN','annlen',FileDescID);
	hdfwarn(len)
	[desc{i},status] = hdf('AN','readann',FileDescID,len);
	hdfwarn(status)
	status = hdf('AN','endaccess',FileDescID);
	hdfwarn(status) 
      end
    end
  end
end

%================================================================
function used = isused(tag,ref,used)
if isempty(used.Tag) && isempty(used.Ref)
  used = 0;
else
  tagIdx = find(used.Tag == tag );
  if isempty(tagIdx)
    used = 0;
  else
    refIdx = find(used.Ref == ref);
    if isempty(refIdx)
      used = 0;
    else
      used = 1;
    end
  end
end

%================================================================
function hdfwarn(status)
%Use the last error generated by the HDF library as a warning.
if status==-1
	msg = '';
	error_code = hdf('HE','value',1);
	msg = sprintf('The NCSA HDF library reported the following error: \n%s',hdf('HE','string',error_code));
	warning(msg)
end

%================================================================
function vdinfo = hdfvdatainfo(filename,fileID,anID,dataset)
%VDATAINFO Information about Vdata set
%
%   VDINFO = VDATAINFO(FILENAME,FILEID,DATASET) returns a structure whose fields
%   contain information about an HDF Vdata data set.  FILENAME is a string
%   that specifies the name of the HDF file.  FILEID is the file identifier
%   returned by using fileid = hdfh('open',filename,permission). DATASET is
%   a string specifying the name of the Vdata or a number specifying the
%   zero based reference number of the data set. If DATASET is the name of
%   the data set and multiple data sets with that name exist in the file,
%   the first dataset is used.  
%
%   Assumptions: 
%               1.  The file has been open.  FILEID is a valid file
%                   identifier.
%               2.  The V and AN interfaces have been started.
%               3.  anID may be -1
%               4.  The file, V, and AN interfaces will be closed elsewhere.
%
%   The fields of VDINFO are:
%
%   Filename          A string containing the name of the file
%   Name              A string containing the name of the data set
%   DataAttributes    An array of structures with fields 'Name' and 'Value'
%                     describing the name and value of the attributes of the
%                     entire data set
%   Class             A string containing the class name of the data set
%   Fields            An array of structures with fields 'Name' and
%                     'Attributes' describing the fields of the Vdata
%   NumRecords        A number specifying the number of records of the data set   
%   IsAttribute       1 if the Vdata is an attribute, 0 otherwise
%   Label             A cell array containing an Annotation label
%   Description       A cell array containing an Annotation description
%   Type              A string describing the type of HDF object 
%   Ref               The reference number of the Vdata set

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.4 $  $Date: 2002/06/05 20:09:08 $

vdinfo = [];

if ~ischar(filename)
  error('FILENAME must be a string.');
end

msg ='Problem connecting to Vdata set. Data set may not exist or file may be corrupt.';
if isnumeric(dataset)
  ref = dataset;
elseif ischar(dataset)
  ref = hdf('VS','find',fileID,dataset);
  if ref == 0
    warning(msg)
    return
  end
else
  error('DATASET must be a number or a string.');
end

vdID = hdf('VS','attach',fileID,ref,'r');
if vdID == -1
  warning(msg);
  return
end

[class, status] = hdf('VS','getclass',vdID);
hdfwarn(status) 

%Get Vdata information
[records, interlace, fieldListLong, size, name, status] = hdf('VS','inquire',vdID);
hdfwarn(status) 
if status==-1
  name = '';
  records = '';
  fields = '';
else
  %Parse field names
  fieldnames = parselist(fieldListLong);
  fields =  cell2struct(fieldnames,'Name',1);

  attrcount = 0;
  for i = 1:length(fieldnames)  
    nfattrs = hdf('VS','fnattrs',vdID,i);
    hdfwarn(nfattrs)
    for j=1:nfattrs
      [name,dataType,count,size,status] = hdf('VS','attrinfo',vdID,i,j);
      hdfwarn(status)
      [attrdata, status] = hdf('VS','getattr',vdID,i,j);
      hdfwarn(status) 
      fields(i).Attributes(j).Name = name;
      fields(i).Attributes(j).Value = attrdata;
      attrcount = attrcount+1;
    end
    if nfattrs==0
      fields(i).Attributes = [];
    end
  end
end

%Get general Vdata attributes
nattrs = hdf('VS','nattrs',vdID);
hdfwarn(nattrs) 
if nattrs>0
  for i = 1:(nattrs-attrcount)
    [name,dataType,count,size,status] = hdf('VS','attrinfo',vdID,'vdata',i);  
    hdfwarn(status) 
    [attrdata, status] = hdf('VS','getattr',vdID,'vdata',i);
    hdfwarn(status) 
    DataAttributes(i).Name = name;
    DataAttributes(i).Value = attrdata;
  end
else
  DataAttributes = [];
end
isattr = logical(hdf('VS','isattr',vdID));

%Get annotations
tag = hdf('ML','tagnum','DFTAG_VS');
if anID ~= -1
  [label,desc] = hdfannotationinfo(filename,anID,tag,ref);
end

%Detach from data set
status = hdf('VS','detach',vdID);
hdfwarn(status)

%Populate output structure
vdinfo.Filename = filename;
vdinfo.Name = name;
vdinfo.Class = class;
vdinfo.Fields = fields;
vdinfo.NumRecords = records;
vdinfo.IsAttribute = isattr;
vdinfo.DataAttributes = DataAttributes;
vdinfo.Label = label;
vdinfo.Description = desc;
%  vdinfo.Tag = hdfml('tagnum','DFTAG_VS');
vdinfo.Ref = ref;
vdinfo.Type = 'Vdata set';

%================================================================
function rinfo = hdfraster8info(filename,imid,anID)
%HDFRASTER8INFO Information about HDF 8-bit Raster image
%
%   RINFO=RASTER8INFO(FILENAME,IMID) returns a structure whose fields contain
%   information about an 8-bit raster image in an HDF file.  FILENAME
%   is a string that specifies the name of the HDF file.  IMID is a string
%   specifying the name of the raster image or a number specifying the
%   image's reference number.  
%
%   The fields of RINFO are:
%
%   Filename       A string containing the name of the file
%   Name           A string containing the name of the image
%   Width          An integer indicating the width of the image in pixels
%   Height         An integer indicating the height of the image in pixels
%   HasPalette     1 if the image has an associated palette, 0 otherwise
%   Ref            Reference number of the raster image
%   Label          A cell array containing an Annotation label
%   Description    A cell array containing an Annotation description
%   Type           A string describing the type of HDF object 

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.3 $  $Date: 2002/06/05 20:09:02 $

rinfo = [];

%Check input arguments
error(nargchk(3,3,nargin));

if ~ischar(filename),	error('FILENAME must be a string.'),	end
if ~hdf('H','ishdf',filename),	error('Invalid HDF file.'),		end
if ~isnumeric(imid),	error('REF must be a number.'),			end

% Chose RIG because the annotations seemed to be linked to this tag.  The
% raster images could be described with all tags, even obsolete ones.
tag = hdf('ML','tagnum','DFTAG_RIG');

status = hdf('DFR8','readref',filename,imid);
if status == -1
  warning('Unable to read image.  The image may not exist or the file may be corrupt.');
else
  [width, height, hasMap, status] = hdf('DFR8','getdims',filename);
  hdfwarn(status);

  %Get annotations
  [label,desc] = hdfannotationinfo(filename,anID,tag,imid);
  
  % Use reference number to name the image: "8-bit Raster Image #refnum".
  % Other browsers use the first data label as the name if it exists. 
  
  name = ['8-bit Raster Image #' num2str(imid)];
  
  %Populate output structure
  rinfo.Filename = filename;
  rinfo.Name = name;
  rinfo.Ref = imid;
  rinfo.Width = width;
  rinfo.Height = height;
  rinfo.HasPalette = hasMap;
  rinfo.Label = label;
  rinfo.Description = desc;
  rinfo.Type = '8-Bit Raster Image';
end

%================================================================
function raster = hdfraster24info(filename, imid, anID)
%HDFRASTER8INFO Information about HDF 24-bit Raster image
%
%   RINFO=RASTER8INFO(FILENAME,IMID) returns a structure whose fields contain
%   information about an 24-bit raster image in an HDF file.  FILENAME
%   is a string that specifies the name of the HDF file.  IMID is a string
%   specifying the name of the raster image or a number specifying the
%   image's reference number.  
%
%   The fields of RINFO are:
%
%   Filename       A string containing the name of the file
%   Name           A string containing the name of the image
%   Width          An integer indicating the width of the image in pixels
%   Height         An integer indicating the height of the image in pixels
%   Interlace      A string describing the interlace mode of the image
%   Ref            Reference number of the raster image
%   Label          A cell array containing an Annotation label
%   Description    A cell array containing an Annotation description
%   Type           A string describing the type of HDF object 

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.3 $  $Date: 2002/06/05 20:09:00 $

rinfo = [];

%Check input arguments
error(nargchk(3,3,nargin));

if ~ischar(filename),	error('FILENAME must be a string.'),	end
if ~hdf('H','ishdf',filename),	error('Invalid HDF file.'),		end
if ~isnumeric(imid),	error('IMID must be a number.'),		end


%Get image information
tag =hdf('ML','tagnum','DFTAG_RI');
status = hdf('DF24','readref',filename,imid);
if status == -1
  warning('Unable to read image.  The image may not exist or the file may be corrupt.');
else
  [width, height, interlace, status] = hdf('DF24','getdims',filename);
  hdfwarn(status);

  [label, desc] = hdfannotationinfo(filename,anID,tag,imid);
  
  % Use reference number to name the image: "24-bit Raster Image #refnum".
  % Other browsers use the first data label as the name if it exists. 
    
  name = ['24-bit Raster Image #' num2str(imid)];

  %Populate output structure
  raster.Filename = filename;
  raster.Name = name;
  raster.Tag = tag;
  raster.Ref = imid;
  raster.Width = width;
  raster.Height = height;
  raster.Interlace = interlace;
  raster.Label = label;
  raster.Description = desc;
  raster.Type = '24-Bit Raster Image';
end

%================================================================
function list = parselist(listin)
%Parse comma seperated list into a cell array
i=0;
while ~isempty(listin)
  i=i+1;
  [list{i}, listin] = strtok(listin,',');
end

%================================================================
function gridinfo = hdfgridinfo(filename,fileID,gridname)
%HDFGRIDINFO Information about HDF-EOS Grid data.
%
%   GRIDINFO = HDFGRIDINFO(FILENAME,GRIDNAME) returns a structure whose
%   fields contain information about a Grid data set in an HDF-EOS
%   file. FILENAME is a string that specifies the name of the HDF-EOS file
%   and GRIDNAME is a string that specifyies the name of the Grid data set
%   in the file.
%
%   The fields of GRIDINFO are:
%
%   Filename       A string containing the name of the file
%   Name           A string containing the name of the Grid
%   UpperLeft      A number specifying the upper left corner location in meters
%   LowerRight     A number specifying the lower right corner location in meters
%   Rows           An integer specifying the number of rows in the Grid
%   Columns        An integer specifying the number of columns in the Grid
%   DataFields     An array of structures with fields 'Name', 'Rank', 'Dims',
%                  'NumberType', 'FillValue', and 'TileDims'. Each structure
%                  describes a data field in the Grid 
%   Attributes     An array of structures with fields 'Name' and 'Value'
%                  describing the name and value of the attributes of the Grid
%   Projection     A structure with fields 'ProjCode', 'ZoneCode',
%                  'SphereCode', and 'ProjParam' describing the Projection
%                  Code, Zone Code, Sphere Code and projection parameters of the Grid
%   Origin Code    A number specifying the origin code for the Grid
%   PixRegCode     A number specifying the pixel registration code
%   Type           A string describing the type of HDF/HDF-EOS
%                  object. 'HDF-EOS Grid' for Grid data sets

%Return empty for data set not found 
gridinfo = [];

error(nargchk(3,3,nargin));

if ~ischar(filename)
  error('FILENAME must be a string.');
end

if ~ischar(gridname)
  error('GRIDNAME must be a string.');
end

%Open interfaces, return early if opening the file or attaching to the grid
%fails
msg = sprintf('Unable to open Grid interface for ''%s'' data set. File may be corrupt.',gridname);
gridID = hdf('GD','attach',fileID,gridname);
if gridID==-1
  warning(msg);
  return
end

%Get upper left and lower right grid corners, # of rows and cols 
[Columns, Rows, UpperLeft, LowerRight, status] = hdf('GD','gridinfo',gridID);
hdfwarn(status) 

%Get info on data fields. fieldListLong is a comma seperated list and
%numberType is a cell array of strings
[nfields, fieldListLong, ranktmp, numberType] = hdf('GD','inqfields',gridID);
hdfwarn(nfields)
if nfields>0
  fieldList = parselist(fieldListLong);
  for i=1:nfields
    fill = hdf('GD','getfillvalue',gridID,fieldList{i});
    [rank{i},dimSizes,numberType{i},dimListLong,status] = hdf('GD','fieldinfo',gridID,fieldList{i});
    dimList = parselist(dimListLong);
    Dims{i} = struct('Name',dimList(:),'Size',num2cell(dimSizes(:)));
    FillValue{i} = fill;
    %Get tile info 
    [tilecode{i},tiledims{i},tilerank,status] = hdf('GD','tileinfo',gridID,fieldList{i});
    if isempty(tiledims{i}) 
      tiledims{i} = [];
    end
  end
  DataFields = struct('Name',fieldList(:),'Rank',rank(:),'Dims',Dims(:),...
		      'NumberType',numberType(:), 'FillValue', FillValue(:),...
		      'TileDims',tiledims(:));
else 
  DataFields = [];
end

%Get attribute information 
[nattrs, attrListLong] = hdf('GD','inqattrs',gridID);
hdfwarn(nattrs)
if nattrs>0
  attrList = parselist(attrListLong);
  Attributes = cell2struct(attrList,'Name',1);
  for i=1:nattrs
    [Attributes(i).Value, status] = hdf('GD','readattr',gridID,attrList{i});
    hdfwarn(status)
  end
else
  Attributes = [];
end

%Get projection information
[Projection.ProjCode,Projection.ZoneCode,Projection.SphereCode,Projection.ProjParam, status] = hdf('GD','projinfo',gridID);
hdfwarn(status)

%Get origin code
OriginCode = hdf('GD','origininfo',gridID);
hdfwarn(OriginCode)

%Get pixel region info
PixRegCode = hdf('GD','pixreginfo',gridID);
hdfwarn(PixRegCode)

%Close interfaces
status = hdf('GD','detach',gridID);

%Assign output structure
gridinfo.Filename     =	 filename;	 
gridinfo.Name         =	 gridname;	 
gridinfo.UpperLeft    =	 UpperLeft;	 
gridinfo.LowerRight   =	 LowerRight;	 
gridinfo.Rows         =	 Rows;	 
gridinfo.Columns      =	 Columns;	 
gridinfo.DataFields   =	 DataFields;	 
gridinfo.Attributes   =	 Attributes;	 
gridinfo.Projection   =	 Projection;	 
gridinfo.OriginCode   =  OriginCode;	 
gridinfo.PixRegCode   =	 PixRegCode;	 
gridinfo.Type         =  'HDF-EOS Grid';

%================================================================
function pointinfo = hdfpointinfo(filename,pointFID,pointname)
%HDFPOINTINFO Information about HDF-EOS Point data.
%
%   POINTINFO = HDFPOINTINFO(FILENAME,POINTNAME) returns a structure whose
%   fields contain information about a Point data set in an HDF-EOS
%   file. FILENAME is a string that specifies the name of the HDF-EOS file
%   and POINTNAME is a string that specifyies the name of the Point data set.
%
%   The fields of POINTINFO are:
%
%   Filename       A string containing the name of the file
%   Name           A string containing the name of the data set
%   Level          An array of structures with fields 'Name', 'NumRecords',
%                  'FieldNames', 'DataType' and 'Index'.  Each structure
%                  describes a level of the Point
%   Attributes     An array of structures with fields 'Name' and 'Value'
%                  describing the name and value of the attributes of the Point
%   Type           A string describing the type of HDF/HDF-EOS object 

%   Assumptions:
%               1.  File has been opened using PT interface.
%               2.  PT interface will be closed elsewhere.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.3 $  $Date: 2002/06/05 20:08:57 $

pointinfo = [];

error(nargchk(3,3,nargin));
if ~ischar(filename)
  error('FILENAME must be a string.');
end
if ~ischar(pointname)
  error('POINTNAME must be a string.');
end

%Open Point interface 
msg = sprintf('Problem retrieving information for Point ''%s''.  File may be corrupt.', pointname);

%Attach to Point
pointID = hdf('PT','attach',pointFID,pointname);
if pointID==-1
  warning(msg);
  return
end

%Get number of levels
NumLevels = hdf('PT','nlevels',pointID);
hdfwarn(NumLevels)

if NumLevels>0
  for i=1:NumLevels
    %Get level name
    [Level(i).Name, status] = hdf('PT','getlevelname',pointID,i-1);
    hdfwarn(status)
    if status==0
      Level(i).Index = hdf('PT','levelindx',pointID,Level(i).Name);
      hdfwarn(Level(i).Index)
    else
      Level(i).Index = [];
    end
    
    %Get number of records
    Level(i).NumRecords = hdf('PT','nrecs',pointID,i-1);
    hdfwarn(Level(i).NumRecords)
    [numfields,fieldsLong,Level(i).DataType,fieldorder] = hdf('PT','levelinfo',pointID,i-1);
    hdfwarn(numfields)
    
    if numfields>0
      fields = parselist(fieldsLong);
      Level(i).FieldNames = fields;
    else 
      Level(i).FieldNames = cell(0);
    end
  end
else
  Level = [];
end

%Get attribute information
[nattrs, attrListLong] = hdf('PT','inqattrs',pointID);
hdfwarn(nattrs)

if nattrs>0
  attrList = parselist(attrListLong);
  Attributes = cell2struct(attrList,'Name',1);
  for i=1:nattrs
    [Attributes(i).Value,status] = hdf('PT','readattr',pointID,attrList{i});
    hdfwarn(status)
  end
else
  Attributes = [];
end

%Close interfaces
status = hdf('PT','detach',pointID);

% Populate output structure
pointinfo.Filename = filename;
pointinfo.Name = pointname;
pointinfo.Attributes = Attributes;
pointinfo.Level = Level;
pointinfo.Attributes = Attributes;
pointinfo.Type = 'HDF-EOS Point';

%================================================================
function swathinfo = hdfswathinfo(filename,fileID,swathname)
%HDFSWATHINFO Information about HDF-EOS Swath data.
%
%   SWATHINFO = HDFSWATHINFO(FILENAME,SWATHNAME) returns a structure whose
%   fields contain information about a Swath data set in an HDF-EOS
%   file. FILENAME is a string that specifies the name of the HDF-EOS file
%   and SWATHNAME is a string that specifies the name of the Swath data set.
%
%   The fields of SWATHINFO are:
%
%   Filename           A string containing the name of the file
%   Name               A string containing the name of the data set
%   DataFields         An array of structures with fields 'Name', 'Rank', 'Dims',
%                      'NumberType', and 'FillValue'.  Each structure
%                      describes a Data field in the Swath 
%   GeolocationFields  An array of structures with fields 'Name', 'Rank', 'Dims',
%                      'NumberType', and 'FillValue'.  Each structure
%                      describes a Geolocation field in the Swath 
%   MapInfo            A structure with fields 'Map', 'Offset', and
%                      'Increment' describing the relationship between the
%                      data and geolocation fields. 
%   IdxMapInfo         A structure with 'Map' and 'Size' describing the
%                      relationship between the indexed elements of the
%                      geolocation mapping
%   Attributes         An array of structures with fields 'Name' and 'Value'
%                      describing the name and value of the swath attributes 
%   Type               A string describing the type of HDF/HDF-EOS object.
%                     'HDF-EOS Swath' for Swath data sets

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.4 $  $Date: 2002/06/05 20:09:06 $

swathinfo = [];

error(nargchk(3,3,nargin));
if ~ischar(filename)
  error('FILENAME must be a string.');
end
if ~ischar(swathname)
  error('SWATHNAME must be a string.');
end

msg = sprintf('Error opening Swath interface to ''%s''.  The Swath may not exist or the file may be corrupt.',swathname);
%Open Swath interfaces
swathID = hdf('SW','attach',fileID,swathname);
if swathID==-1
  warning(msg);
  return
end

%Get Data Field information
[nfields, fieldListLong, ranktmp, numbertype] = hdf('SW','inqdatafields',swathID);
hdfwarn(nfields)
if nfields>0
  fieldList = parselist(fieldListLong);
  for i=1:nfields
    [fill, status] = hdf('SW','getfillvalue',swathID,fieldList{i});
%    hdfwarn(status);
    [rank{i},dimSizes,numberType{i},dimListLong,status] = hdf('SW','fieldinfo',swathID,fieldList{i});
    dimList = parselist(dimListLong);
    Dims{i} = struct('Name',dimList(:),'Size',num2cell(dimSizes(:)));
    FillValue{i} = fill;
  end
  DataFields = struct('Name',fieldList(:),'Rank',rank(:),'Dims',Dims(:),...
		      'NumberType',numberType(:),'FillValue', FillValue(:));
else
  DataFields = [];
end

%Get Geolocation information
Dims = {};
FillValue = {};
rank = {};

[ngeofields, fieldListLong, ranktmp, numbertype] = hdf('SW','inqgeofields',swathID);
hdfwarn(ngeofields)
if ngeofields>0
  fieldList = parselist(fieldListLong);
  for i=1:ngeofields
    [fill, status] = hdf('SW','getfillvalue',swathID,fieldList{i});
    %    hdfwarn(status);
    [rank{i},dimSizes,numberType,dimListLong,status] = hdf('SW','fieldinfo',swathID,fieldList{i});
    dimList = parselist(dimListLong);
    Dims{i} = struct('Name',dimList(:),'Size',num2cell(dimSizes(:)));
    FillValue{i} = fill;
  end
  GeolocationFields = struct('Name',fieldList(:),'Rank',rank(:),'Dims',Dims(:),'NumberType',numbertype(:),'FillValue',FillValue(:));
else
  GeolocationFields = [];
end

%Get Geolocation relations
[nmaps, dimMapLong, offset, increment] = hdf('SW','inqmaps',swathID);
if nmaps>0
  dimMap = parselist(dimMapLong);
  MapInfo = struct('Map',dimMap(:),'Offset',num2cell(offset(:)),'Increment',num2cell(increment(:)));
else 
  MapInfo = [];
end

%Get index mapping relations
[nmaps, idxMapLong, idxSizes] = hdf('SW','inqidxmaps',swathID);
hdfwarn(nmaps)
if nmaps>0
  idxMap = parselist(idxMapLong);
  IdxMapInfo = struct('Map',idxMap(:),'Size',num2cell(idxSizes(:)));
else 
  IdxMapInfo = [];
end

%Retrieve attribute information
[nattr, attrListLong] = hdf('SW','inqattrs',swathID);
hdfwarn(nattr)
if nattr>0
  attrList = parselist(attrListLong);
  Attributes = cell2struct(attrList,'Name',1);
  for i=1:nattr
    [Attributes(i).Value, status] = hdf('SW','readattr',swathID,attrList{i});
    hdfwarn(status)
  end
else 
  Attributes = [];
end

%Close interfaces
status = hdf('SW','detach',swathID);

%Populate output structure
swathinfo.Filename         = filename;
swathinfo.Name             = swathname;         
swathinfo.DataFields       = DataFields;       
swathinfo.GeolocationFields= GeolocationFields;
swathinfo.MapInfo          = MapInfo;
swathinfo.IdxMapInfo       = IdxMapInfo;
swathinfo.Attributes       = Attributes;
swathinfo.Type             = 'HDF-EOS Swath';

%================================================================
function [vginfo, children] = hdfvgroupinfo(filename,dataset,fileID,sdID,anID)
%HDFVGROUPINFO Information about HDF Vgroup
%
%   VGINFO = HDFVGROUPINFO(FILENAME,DATASET) returns a structure whose
%   fields contain information about an HDF Vgroup.  FILENAME is a string
%   that specifies the name of the HDF file.   DATASET is a string
%   specifying the name of the Vgroup (?) or a number specifying the zero based
%   index of the data set.  FILEID is the file identifier returned by
%   hdfh('open',... and SDID is the sds identifier returned by
%   hdfsd('start',...
%
%   Assumptions: 
%               1.  The file has been open.  FILEID is a valid file
%                   identifier.
%               2.  The V, AN and SD interfaces have been started.
%               3.  anID may be -1.
%               4.  The V, AN and SD interface and file will be closed elsewhere.
% 
%   The fields of VGINFO depend on the contents of the file.
%
%   Possible fields of VGINFO are:
%
%   Filename  A string containing the name of the file
%   Name      A string containing the name of the data set
%   Class     A string containing the name of the class of the data set
%   Vgroup    An array of structures describing Vgroups
%   SDS       An array of structures describing Scientific Data sets
%   Vdata     An array of structures describing Vdata sets
%   Raster24  An array of structures describing 24-bit raster images  
%   Raster8   An array of structures describing 8-bit raster images
%   Type      A string describing the type of HDF object 
%   Tag       The tag number of this Vgroup
%   Ref       The reference number of this vgroup
%

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.4 $  $Date: 2002/06/05 20:09:10 $

%Check input arguments
error(nargchk(5,5,nargin));

if (~ischar(filename)),		error('FILENAME must be a string.'),	end

%if ~isnumeric(dataset)
%  error('DATASET must be a number.');
%end

%Initialization 
vginfo = struct('Filename',[],'Name',[],'Class',[],'Vgroup',[],'SDS',[],'Vdata',[],'Raster24',[],'Raster8',[],'Tag',[],'Ref',[]);
vgroup = [];
sds = [];
raster8 = [];
raster24 = [];
vdata = [];
children.Tag = [];
children.Ref = [];

msg ='Problem connecting to Vgroup. File may be corrupt or Vgroup may not exist.';
vgID = hdf('V','attach',fileID,dataset,'r');
if vgID == -1
  warning(msg);
  return
end

%Get name and class of Vgroup
[numEntries, name, status] = hdf('V','inquire',vgID);
hdfwarn(status)
[class, status] = hdf('V','getclass',vgID);
hdfwarn(status)

%Get child/tag reference pairs
childTags = [];
childRefs = [];
numChildren = hdf('V','ntagrefs',vgID);
hdfwarn(numChildren)
if numChildren>0
  [childTags, childRefs, status] = hdf('V','gettagrefs',vgID,numChildren);
  hdfwarn(status)
  children.Tag = childTags;
  children.Ref = childRefs;
end

%Close interfaces
status = hdf('V','detach',vgID);
hdfwarn(status)

%Get children information. Classes of Vgroups to ignore: cdf0.0, rig0.0
if ~strcmp(lower(class),'cdf0.0') & ~strcmp(lower(class),'rig0.0')
  %Use same tag mapping as NCSA
  for i=1:length(childTags)
    switch childTags(i)
     case {hdf('ML','tagnum','DFTAG_SD'),hdf('ML','tagnum','DFTAG_NDG'),hdf('ML','tagnum','DFTAG_SDG')}
      %Change reference number to index
      index = hdf('SD','reftoindex',sdID,childRefs(i));
      hdfwarn(index)
      if index~=-1
	sds = [sds hdfsdsinfo(filename,sdID,anID,index)]; 
      end
     case hdf('ML','tagnum','DFTAG_RI8')
      raster8 = [raster8 hdfraster8info(filename,anID,childRefs(i))];
     case hdf('ML','tagnum','DFTAG_RI')
      raster24 = [raster24 hdfraster24info(filename,anID,childRefs(i))];
     case hdf('ML','tagnum','DFTAG_RIG')
      raster24 = [raster24 hdfraster24info(filename,anID,childRefs(i))];
     case {hdf('ML','tagnum','DFTAG_VS'),hdf('ML','tagnum','DFTAG_VH')}
      vdata = [vdata hdfvdatainfo(filename,fileID,anID,childRefs(i))];
     case hdf('ML','tagnum','DFTAG_VG')
      [vgrouptemp, child] = hdfvgroupinfo(filename,childRefs(i),fileID,sdID,anID);
      vgroup = [vgroup vgrouptemp];
      if ~isempty(child.Tag)
	children.Tag = [children.Tag child.Tag];
	children.Ref = [children.Ref child.Ref];
      end
     otherwise
      warning(['Unknown child in Vgroup with Tag number: ' num2str(childTags(i))]);
    end
  end               
  %Populate output structure
  vginfo.Filename = filename;
  vginfo.Name = name;
  vginfo.Class = class;
  if ~isempty(vgroup),		vginfo.Vgroup = vgroup;		end
  if ~isempty(sds),			vginfo.SDS = sds;			end
  if ~isempty(vdata),		vginfo.Vdata = vdata;		end
  if ~isempty(raster24),	vginfo.Raster24 = raster24;	end
  if ~isempty(raster8),		vginfo.Raster8 = raster8;	end
  vginfo.Tag = hdf('ML','tagnum','DFTAG_VG');
  vginfo.Ref = dataset;
  vginfo.Type = 'Vgroup';
end

%================================================================
function [sdinfo, IsScale] = hdfsdsinfo(filename, sdID, anID, dataset)
%HDFSDSINFO Information about HDF Scientific Data Set. 
%
%   [SDINFO, IsScale]=SDSINFO(FILENAME,DATASET) returns a structure SDINFO whose
%   fields contain information about an HDF Scientific Data
%   Set(SDS). IsScale is a true (1) if the SDS is a dimension scale, false
%   (0) otherwise. FILENAME is a string that specifies the name of the HDF
%   file.  SDID is the sds identifier returned by hdfsd('start',... ANID is
%   the annotation identifier returned by hdfan('start',... DATASET is a
%   string specifying the name of the SDS or a number specifying the zero
%   based index of the data set. If DATASET is the name of the data set and
%   multiple data sets with that name exist in the file, the first dataset
%   is used.
%
%   Assumptions: 
%               1.  The file has been open.
%               2.  The SD and AN interfaces have been started.  
%               3.  anID may be -1
%               3.  The SD and AN interfaces and file will be closed elsewhere.
%
%   The fields of SDINFO are:
%
%   Filename          A string containing the name of the file
%   Type              A string describing the type of HDF object 
%   Name              A string containing the name of the data set
%   Rank              A number specifying the number of dimensions of the data set
%   DataType          A string specifying the precision of the data
%   Attributes        An array of structures with fields 'Name' and 'Value'
%                     describing the name and value of the attributes of the data set
%      Dims           An array of structures with fields 'Name', 'DataType',
%                     'Size', 'Scale', and 'Attributes'.  Describing the the
%                     dimensions of the data set.  'Scale' is an array of numbers
%                     to place along the dimension and demarcate intervals in
%                     the data set.
%   Label             A cell array containing an Annotation label
%   Description       A cell array containing an Annotation descriptoin
%   Index             Number indicating the index of the SDS

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.4 $  $Date: 2002/06/05 20:09:04 $

sdinfo = struct('Filename',[],'Type',[],'Name',[],'Rank',[],'DataType',[],'Attributes',[],'Dims',[],'Label',[],'Description',[],'Index',[]);

error(nargchk(4,4,nargin));

if ~ischar(filename)
  error('FILENAME must be a string.');
end

msg ='Problem connecting to Scientific Data set. File may be corrupt or data set may not exist.';
%User may input name or index to the SDS
if isnumeric(dataset)
  index = dataset;
elseif ischar(dataset)
  index = hdf('SD','nametoindex',sdID,dataset);
  if index == -1
    warning(msg);
    return
  end
else
  error('DATASET must be a number or a string.');
end

sdsID = hdf('SD','select',sdID,index);
if sdsID == -1
  warning(msg);
  return
end

%Convert index to reference number 
ref = hdf('SD','idtoref',sdsID);
hdfwarn(ref)

%Get lots of info
[sdsName, rank, dimSizes, sddataType, nattrs, status] = hdf('SD','getinfo',sdsID);
hdfwarn(status)

%Get SD attribute information. The index for readattr is zero based.
if nattrs>0
  for i = 1:nattrs
    [arrayAttribute(i).Name,attrDataType,count,status] = hdf('SD','attrinfo',sdsID,i-1);
    hdfwarn(status)
    [arrayAttribute(i).Value, status] = hdf('SD','readattr',sdsID,i-1);
    hdfwarn(status)
  end
else
  arrayAttribute = [];
  sdDataType = '';
end
  
IsScale = logical(hdf('SD','iscoordvar',sdsID));

%If it is not a dimension scale, get dimension information
%Dimension numbers are 0 based (?)
if IsScale == 0
  for i=1:rank
    dimID = hdf('SD','getdimid',sdsID,i-1);
    %Use sizes from SDgetinfo because this size may be Inf
    [dimName{i}, sizeDim,DataType{i}, nattrs, status] = hdf('SD','diminfo',dimID);
    hdfwarn(status)
    if strcmp(DataType{i},'none')
      Scale{i} = 'none';
    elseif isinf(sizeDim)
      Scale{i} = 'unknown';
    else
      try
        [Scale{i},status] = hdf('SD','getdimscale',dimID);
      catch
        Scale{i} = 'none';
      end
    end
    Size{i} = dimSizes(i);
    if nattrs>0
      for j=1:nattrs
	[Name{j},dataType,count,status] = hdf('SD','attrinfo',dimID,j-1);
	hdfwarn(status)
	[Value{j}, status] = hdf('SD','readattr',dimID,j-1);
	hdfwarn(status)
      end
      Attributes{i} = struct('Name',Name(:),'Value',Value(:));
    else
      Attributes = [];
    end
  end
  dims = struct('Name',dimName(:),'DataType',DataType(:),'Size',Size(:),'Scale',Scale(:),'Attributes',Attributes(:));
else
  dims = [];
end

%Get any associtated annotations
tag = hdf('ML','tagnum','DFTAG_NDG');
if anID ~= -1
  [label,desc] = hdfannotationinfo(filename,anID,tag,ref);
  if isempty(label) || isempty(desc)
    tag = hdf('ML','tagnum','DFTAG_SD');
    [label,desc] = hdfannotationinfo(filename,anID,tag,ref);
  end
end

%Close interfaces
status = hdf('SD','endaccess',sdsID);
hdfwarn(status)

%Populate output structure
sdinfo.Filename = filename;
sdinfo.Name = sdsName;
sdinfo.Index = index;
sdinfo.Rank = rank;
sdinfo.DataType = sddataType;
if ~isempty(arrayAttribute)
  sdinfo.Attributes = arrayAttribute;
end
if ~isempty(dims)
  sdinfo.Dims = dims;
end
sdinfo.Label = label;
sdinfo.Description = desc;
%sdinfo.IsScale = IsScale;
sdinfo.Type = 'Scientific Data Set';

%================================================================
function [label, desc] = hdfannotationinfo(filename,anID,tag,ref)
%HDFANNOTATIONINFO Retrieve about HDF Annotation. 
%
%   [LABEL, DESC] = HDFANNOTATIONINFO(FILENAME,TAG,REF) returns the data
%   description DESC, and data label LABEL, for the HDF object described by
%   the TAG, REF pair in the file FILENAME.  FILEID is the file identifier
%   returned by hdfh('open',... and ANID is the AN identifier returned by
%   hdfan('start',...  If no label or description exist for the HDF object,
%   then LABEL and DESC will be empty cell arrays.
%
%   Assumptions: 
%               1.  The file has been open.  FILEID is a valid file
%                   identifier.
%               2.  The AN interface has been started.
%               3.  The AN interface and file will be closed elsewhere.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.3 $  $Date: 2002/06/05 20:08:53 $

label = {};		desc = {};

% hdfan('annlen',...) does not include null termination character read by 
% hdfan('readann'...). Adding 1 to length for this reason.
numDataLabel = hdf('AN','numann',anID,'data_label',tag,ref);
hdfwarn(numDataLabel)
if numDataLabel>0
  [DataLabelID,status] = hdf('AN','annlist',anID,'data_label',tag,ref);
  hdfwarn(status)
  if status~=-1
    for i=1:numDataLabel
      length = hdf('AN','annlen',DataLabelID(i));
      hdfwarn(length)
      [label{i},status] = hdf('AN','readann',DataLabelID(i),length+1);
      hdfwarn(status)
    end
  end
end
numDataDesc = hdf('AN','numann',anID,'data_desc',hdf('ML','tagnum','DFTAG_NDG'),ref);
hdfwarn(numDataDesc)
if numDataDesc >0
  [DataDescID, status] = hdf('AN','annlist',anID,'data_desc',hdf('ML','tagnum','DFTAG_NDG'),ref);
  if status~=-1
    for i=1:numDataDesc
      length = hdf('AN','annlen',DataDescID(i));
      hdfwarn(length)
      desc{i} = hdf('AN','readann',DataDescID(i),length+1);
      hdfwarn(status)
    end
  end
end

% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
function varargout = hdfread(varargin)
%HDFREAD extract data from HDF file
%   
%   HDFREAD reads data from a data set in an HDF or HDF-EOS file.  If the
%   name of the data set is known, then HDFREAD will search the file for the
%   data.  Otherwise, use HDFINFO to obtain a structure describing the
%   contents of the file. The fields of the structure returned by HDFINFO are
%   structures describing the data sets contained in the file.  A structure
%   describing a data set may be extracted and passed directly to HDFREAD.
%   These options are described in detail below.
%   
%   DATA = HDFREAD(FILENAME,DATASETNAME) returns in the variable DATA all 
%   data from the file FILENAME for the data set named DATASETNAME.  
%   
%   DATA = HDFREAD(HINFO) returns in the variable DATA all data from the
%   file for the particular data set described by HINFO.  HINFO is a
%   structure extracted from the output structure of HDFINFO.
%   
%   DATA = HDFREAD(...,PARAMETER,VALUE,PARAMETER2,VALUE2...) subsets the
%   data according to the string PARAMETER which specifies the type of
%   subsetting, and the values VALUE.  The table below outlines the valid
%   subsetting parameters for each type of data set.  Parameters marked as
%   "required" must be used to read data stored in that type of data set.
%   Parameters marked "exclusive" may not be used with any other subsetting
%   parameter, except any required parameters.  When a parameter requires
%   multiple values, the values must be stored in a cell array.  Note that
%   the number of values for a parameter may vary for the type of data set.
%   These differences are mentioned in the description of the parameter.
%
%   [DATA,MAP] = HDFREAD(...) returns the image data and the colormap for an
%   8-bit raster image.
%   
%   Table of available subsetting parameters
%
%
%           Data Set          |   Subsetting Parameters
%          ========================================
%           HDF Data          |
%                             |
%             SDS             |   'Index'
%                             |
%             Vdata           |   'Fields'
%                             |   'NumRecords'
%                             |   'FirstRecord'
%          ___________________|____________________
%           HDF-EOS Data      |   
%                             |
%             Grid            |   'Fields'         (required)
%                             |   'Index'          (exclusive)
%                             |   'Tile'           (exclusive)
%                             |   'Interpolate'    (exclusive)
%                             |   'Pixels'         (exclusive)
%                             |   'Box'
%                             |   'Time'
%                             |   'Vertical'
%                             |
%             Swath           |   'Fields'         (required)
%                             |   'Index'          (exclusive)
%                             |   'Time'           (exclusive)
%                             |   'Box'
%                             |   'Vertical'
%                             |   'ExtMode'
%                             |
%             Point           |   'Level'          (required)
%                             |   'Fields'         (required)
%                             |   'RecordNumbers'
%                             |   'Box'
%                             |   'Time'
%
%    There are no subsetting parameters for Raster Images
%
%
%   Valid parameters and their values are:
%
%   'Index' 
%   Values for 'Index': START, STRIDE, EDGE
%
%     START, STRIDE and EDGE must be arrays the same size as the
%     number of dimensions. START specifies the location in the data set to
%     begin reading.  Each number in START must be smaller than its
%     corresponding dimension.  STRIDE is an array specifying the interval
%     between the values to read.  EDGE is an array specifying the length of
%     each dimension to read.  The region specified by START, STRIDE and EDGE
%     must be within the dimensions of the data set.  If either START, 
%     STRIDE, or EDGE is empty, then default values are calculated assuming:
%     starting at the first element of each dimension, a stride of one, and
%     EDGE to read the from the starting point to the end of the dimension.
%     The defaults are all ones for START and STRIDE, and EDGE is an array
%     containing the lengths of the corresponding dimensions.  START,STRIDE
%     and EDGE are one based. START,STRIDE and EDGE vectors must be stored
%     in a cell as in the following notation: {START,STRIDE,EDGE}.
%
%   'FIELDS'
%    Values for 'Fields' are: FIELDS
%
%      Read data from the field(s) FIELDS of the data set.  FIELDS must be a
%      single string.  For multiple field names, use a comma separated list.
%      For Grid and Swath data sets, only one field may be specified.
%
%   'Box'
%   Values for 'Box' are: LONG, LAT, MODE
%
%     LONG and LAT are numbers specifying a latitude/longitude  region. MODE
%     defines the criterion for the inclusion of a cross track in a region.
%     The cross track in within a region if its midpoint is within the box,
%     either endpoint is within the box or any point is within the box.
%     Therefore MODE can have values of: 'midpoint', 'endpoint', or
%     'anypoint'. MODE is only valid for Swath data sets and will be ignored
%     if specified for Grid or Point data sets.
%
%   'Time'
%   Values for 'Time' are: STARTTIME, STOPTIME, MODE
%
%     STARTTIME and STOPTIME are numbers specifying a region of time. MODE
%     defines the criterion for the inclusion of a cross track in a region.
%     The cross track in within a region if its midpoint is within the box,
%     either endpoint is within the box or any point is within the box.
%     Therefore MODE can have values of: 'midpoint', 'endpoint', or
%     'anypoint'. MODE is only valid for Swath data sets and will be ignored
%     if specified for Grid or Point data sets.
%
%   'Vertical'
%   Values for 'Vertical' are: DIMENSION, RANGE
%
%     RANGE is a vector specifying the min and max range for the
%     subset. DIMENSION is the name of the field or dimension to subset by.  If
%     DIMENSION is the dimension, then the RANGE specifies the range of
%     elements to extract (1 based).  If DIMENSION is the field, then RANGE
%     specifies the range of values to extract. Vertical subsetting may be
%     used in conjunction with 'Box' and/or 'Time'.  To subset a region along
%     multiple dimensions, vertical subsetting may be used up to 8 times in
%     one call to HDFREAD.
%
%   'ExtMode'
%   Values for 'ExtMode' are: EXTMODE
%
%     EXTMODE is either 'Internal' (default) or 'External'.  If the mode is
%     set to 'Internal then the geolocation fields and data fields must be
%     in the same swath.  If the mode is set to 'External' then the
%     geolocation fields and data fields may be in different swaths.  This
%     parameter is only used for Swath data when extracting a time period or
%     a region.
%
%   'Pixels'
%   Values for 'Pixels' are: LON, LAT
%
%     LON and LAT are numbers specifying a latitude/longitude region.  The
%     longitude/latitude region will be converted into pixel rows and
%     columns with the origin in the upper left-hand corner of the grid.
%     This is the pixel equivalent of reading a 'Box' region.
%
%   'RecordNumbers'
%   Available parameter for 'RecordNumbers' is: RecNums
%
%     RecNums is a vector specifying the record numbers to read.  
%
%   'Level'
%   Value for 'Level' is: LVL
%   
%     LVL is a string representing the name of the level to read or a one
%     based number specifying the index of the level to read or  from an
%     HDF-EOS Point data set.
%
%   'NumRecords'
%   Available parameter for 'NumRecords' is: NumRecs
%
%     NumRecs is a number specifying the total number of records to read.
%
%   'FirstRecord'
%   Required value for 'FirstRecord' is: FirstRecord
%
%     FirstRecord is a one based number specifying the first record from which
%     to begin reading.
%
%   'Tile'
%   Required value for 'Tile' is: TileCoords
%
%     TileCoords is a vector specifying the tile coordinates to read. The
%     elements of TileCoords are one based numbers.
%
%   'Interpolate'
%   Values for 'Interpolate' are: LON, LAT
%
%     LON and LAT  are numbers specifying a latitude/longitude
%     points for bilinear interpolation.
%
%    References: 
%
%    Example 1:
%            
%             %  Read data set named 'Example SDS' 
%             data = hdfread('example.hdf','Example SDS');
%
%    Example 2:
%
%             %  Retrieve info about example.hdf
%             fileinfo = hdfinfo('example.hdf');
%             %  Retrieve info about Scientific Data Set in example.hdf
%             data_set_info = fileinfo.SDS;
%             %  Check the size
%             data_set_info.Dims.Size
%             % Read a subset of the data using info structure
%             data = hdfread(data_set_info,...
%                              'Index',{[3 3],[],[10 2 ]});
%
%   See also HDFTOOL, HDFINFO, HDF.  
  
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.7.12.1 $  $Date: 2007/07/24 21:25:10 $

varargout{1} = [];
[hinfo,subsets] = dataSetInfo(varargin{:});

if isempty(hinfo)
  error('hdfread:noDataSets', 'No data set found. Make sure the filename and data set name are correct.');
end

[start,stride,edge,fields,numrecords,firstRecord,level] = parseSubsets(subsets);

switch hinfo.Type
 case 'Scientific Data Set'
  varargout{1} = hdfsdsread(hinfo,start,stride,edge);
 case 'Vgroup'
  error('hdfread:specificDataset', ...
        'A Vgroup is a container for other data sets. A specific data set must be specified.');
 otherwise 
  error('hdfread_m:datatype', ...
        'Data set type not recognized. Verify the HINFO structure is correct.\nConsider using HDFINFO to obtain this structure.');
end

%================================================================
function [start,stride,edge,fields,numrecords,firstRecord,level] = parseSubsets(subsets)
%PARSESUBSETS 
%  Parse some of the subsetting param/value pairs. Values for parameters
%  that are required for data sets are extracted from the variable list of
%  subsetting parameters. This routine will error if the input parameters are
%  not consistent with the param/value syntax described in the help for HDFREAD.

%Return empty structures if not assigned on the command line
start = [];
stride = [];
edge = [];
fields = [];
numrecords = [];
firstRecord = [];
level = [];

if rem(length(subsets),2)
  error('MATLAB:hdfread:subsetValuePairs', 'The subset/value inputs must always occur as pairs.');
end

%Parse subsetting parameters
numPairs = length(subsets)/2;
params = subsets(1:2:end);
values = subsets(2:2:end);

allStrings = {'index','fields','numrecords','firstrecord','tile','interpolate',...
	      'pixels','box','time','vertical','extmode','level','recordnumbers'};
for i=1:length(params)
  idx = strmatch(lower(params{i}),allStrings);
  switch length(idx)
   case 0
    error('MATLAB:hdfread:unknownArgument', 'Unknown string argument: "%s."', params{i});
   case 1
    params{i} = allStrings{idx};
   otherwise
    error('MATLAB:hdfread:ambiguousArgument', 'Ambiguous string argument: "%s."',params{i});
  end
end

cellmsg = '''%s'' method requires %i value(s) to be stored in a cell array.';
for i=1:numPairs
  switch params{i}
   case 'index'
    if iscell(values{i})
      if length(values{i})<3
        error('MATLAB:hdfread:notEnoughInputsToIndex', cellmsg,params{i}, 3);
      else
        [start,stride,edge] = deal(values{i}{:});
      end
    else
      error('MATLAB:hdfread:indexIsNotCell', cellmsg,params{i}, 3);
    end
   case 'fields'
    % 1 comma separated string, 1 cell w/comma separated string, 
    % or 1 cell array of strings are all valid values
    if iscell(values{i})
      if iscellstr(values{i})
		fields = sprintf('%s,',values{i}{:});
		fields = fields(1:end-1);
      end
    else
      fields = values{i};
    end
   case 'numrecords'
    if iscell(values{i})   
      if length(values{i})>1
        error('MATLAB:hdfread:tooManyInputsToNumrecords', cellmsg,params{i}, 1);
      end
      numrecords = values{i}{:};
    else
      numrecords = values{i};
    end
   case 'firstrecord'
    if iscell(values{i})
      if length(values{i})>1
        error('MATLAB:hdfread:tooManyInputsToFirstrecord', cellmsg, params{i}, 1);
      end
      firstRecord = values{i}{:}; 
    else
      firstRecord = values{i};
    end
   case 'level'
    if iscell(values{i})
      if length(values{i})>1
        error('MATLAB:hdfread:tooManyInputsToLevel', cellmsg, params{i}, 1);
      end
      level = values{i}{:};
    else
      level = values{i};
    end
  end
end

%=================================================================
function [hinfo,subsets] = dataSetInfo(varargin)
%DATASETINFO Return info structure for data set and subset param/value pairs
%
%  Distinguish between DATA = HDFREAD(FILENAME,DATASETNAME) and 
%  DATA = HDFREAD(HINFO)

msg = sprintf('Invalid input arguments. HDFREAD requires a filename and data set name, \nor an information structure obtained from HDFINFO.');
if (nargin < 1),		error('MATLAB:hdfread:invalidInput', msg),		end

if ischar(varargin{1}) %HDFREAD(FILENAME,DATASETNAME...)

  msg = nargchk(2,inf,nargin);
  if (~isempty(msg)),	error('MATLAB:hdfread:numberOfInputs', msg);	end
    
  filename = varargin{1};
  %Get full filename
  fid = fopen(filename);
  if fid ~= -1
    filename = fopen(fid);
    fclose(fid);
  else
    error('MATLAB:hdfread:fileOpen', 'File not found.');
  end

  %Use HX interface in case data is in external files
  hdf('HX', 'setdir', fileparts(filename));
  if ischar(varargin{2})
    dataname = varargin{2};
    hinfo = hdfquickinfo(filename,dataname);
    subsets = varargin(3:end);
  elseif isnumeric(varargin{2}) %Obsolete syntax
    subsets = [];
    hinfo.Filename = filename;
    hinfo.TagRef = varargin{2};
    hinfo.Type = 'Obsolete';
    warning('MATLAB:hdfread:obsoleteUsage', ...
            'This usage of HDFREAD is obsolete and may be removed in a future version.\nConsider using IMREAD instead.');
  else
    error('MATLAB:hdfread:invalidDatasetName', sprintf('%s', msg)); %Invalid input
  end
elseif isstruct(varargin{1}) %HDFREAD(HINFO,...)
  hinfo = varargin{1};
  if length(hinfo)>1
    error('MATLAB:hdfread:tooManyHINFO', 'HINFO must be a structure describing a specific data set in the file.');
  end
  if ~isfield(hinfo,'Type')
    error('MATLAB:hdfread:missingTypeInHINFO', ...
          'HINFO is not a valid structure describing an HDF or HDF-EOS data set.  \nConsider using HDFINFO to obtain this structure.');
  end
  subsets = varargin(2:end);
else %Invalid input
  error('MATLAB:hdfread:invalidFilename', sprintf('%s', msg));
end

% ------------------------------------------------------------------------------------
function hinfo = hdfquickinfo(filename,dataname)
%HDFQUICKINFO scan HDF file
%
%   HINFO = HDFQUICKINFO(FILENAME,DATANAME) scans the HDF file FILENAME for
%   the data set named DATANAME.  HINFO is a structure describing the data
%   set.  If no data set is found an empty structure is returned.

%Scientific Data Set
[found, hinfo] = findInsideSD ( filename, dataname );

% ------------------------------------------------------------------------------------
function [found, hinfo] = findInsideSD ( filename, dataname )
% If given something like "/varname", then the slash just means
% that varname is part of the root group.  We need to remove the slash.
if ( dataname(1) == '/' ),		dataname(1) = '';	end
found = false;	hinfo = [];

sdID = hdf('SD','start',filename,'read');
fileID = hdf('H','open',filename,'read',0);
anID = hdf('AN','start',fileID);
index = hdf('SD','nametoindex',sdID,dataname);
if index~=-1
    found = 1;
    hinfo = hdfsdsinfo(filename,sdID,anID,dataname);
end
%Close interface
hdf('SD','end',sdID);
hdf('AN','end',anID);
hdf('H','close',fileID);

% ------------------------------------------------------------------------------------
function [found, hinfo] = findInsideVdata ( filename, dataname )
% If given something like "/varname", then the slash just means
% that varname is part of the root group.  We need to remove the slash.
if ( dataname(1) == '/' ),		dataname(1) = '';	end
found = false;	hinfo = [];

fileID = hdf('H','open',filename,'read',0);
anID = hdf('AN','start',fileID);
hdf('V','start',fileID);
ref = hdf('VS','find',fileID,dataname);
if ref~=0
    found = 1;
    hinfo = hdfvdatainfo(filename,fileID,anID,ref);
end
hdf('V','end',fileID);
hdf('AN','end',anID);
hdf('H','close',fileID);

% ---------------------------------------------------------------------------
function varargout = hdfml(varargin)
%HDFML MATLAB-HDF gateway utilities.
resp = [];
if nargout>0
  [varargout{1:nargout}] = hdf('ML',varargin{:});
else
  hdf('ML',varargin{:})
end
if (~isempty(resp)),	varargout{1} = resp;	end

% ---------------------------------------------------------------------------
function data  = hdfsdsread(hinfo,start,stride,edge)
%HDFSDSREAD read HDF Scientific Data Set%

%Parse inputs and assign default parameters
[start,stride,edge] = parseSDSInputs_hdfsdsread(hinfo,start,stride,edge);

sdID = hdf('SD','start',hinfo.Filename,'read');
if sdID == -1
  error('MATLAB:hdfsdsread:sdsStart', 'Problem opening file %s. The file may be corrupt.',hinfo.File);
end
sdsID = hdf('SD','select',sdID,hinfo.Index);
if sdsID == -1
  hdf('SD','end',sdID);
  error('MATLAB:hdfsdsread:sdsSelect', ...
        'Problem reading Scientific Data Set ''%s''. The data set may not exist or file may be corrupt.',hinfo.Name) ;
end

%  HDFSD('readdata',... will error with incorrect input arguments.  To prevent
%  leaving open identifiers if an error occurs, catch the error then and return
%  a warning.
try
  [data,status] = hdf('SD','readdata',sdsID,start,stride,edge);
catch
  error('MATLAB:HDFSDSREAD:sdsReaddata', lasterr);
end
if status == -1
  error('MATLAB:HDFSDSREEAD:readdata', ...
        'Problem reading Scientific Data Set ''%s''. The data set may not exist or file may be corrupt.',hinfo.Name);
end

%Permute data to be the expected dimensions
data = permute(data,ndims(data):-1:1);

status = hdf('SD','endaccess',sdsID);
if status == -1
    hdf('SD','end',sdID);
    error('MATLAB:HDFSDSREAD:endaccess', 'Could not close dataset %s.', hinfo.Name);
end
status = hdf('SD','end',sdID);
if status == -1
    error('MATLAB:HDFSDSREAD:end', 'Could not close file %s.', hinfo.Filename);
end

%============================================================
function [start,stride,edge] = parseSDSInputs_hdfsdsread(hinfo,start,stride,edge)
%Check for valid inputs to HDFSDSREAD
%There must be START, STRIDE, and EDGE parameters
error(nargchk(1,4,nargin));

msg = 'HINFO is not a valid structure describing a Scientific Data Set.  Consider using HDFINFO to obtain this structure.';
if (~isstruct(hinfo)),	error(msg);		end

%Check for required fields in hinfo structure
fNames = fieldnames(hinfo);
numFields = length(fNames);
reqFields = {'Filename','Rank','Dims','Index','Name'};
numReqFields = length(reqFields);
if numFields >= numReqFields
  for i=1:numReqFields
    if ~isfield(hinfo,reqFields{i})
      error(msg);
    end
  end
else 
  error(msg);
end
if ~isfield(hinfo.Dims,'Size')
  error(msg)
end

%Assign default values to parameters not defined in input
%start, stride and edge are one based. 
if any([start<1, stride<1, edge<1])
  error('START, STRIDE, and EDGE values must be 1 or greater.');
end

if isempty(start)
  start = zeros(1,hinfo.Rank);
else
  start = start-1;
end

if isempty(stride)
  stride = ones(1,hinfo.Rank);
end

if isempty(edge)
  for i=1:hinfo.Rank
    edge(i) = fix((hinfo.Dims(i).Size-start(i))/stride(i));
  end
end
