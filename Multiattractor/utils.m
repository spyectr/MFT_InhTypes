classdef utils
    methods(Static)
        
        
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        %% v2struct
        % v2struct Pack/Unpack Variables to/from a scalar structure.
        function varargout = v2struct(varargin)
            
            %% Description
            %    v2struct has dual functionality in packing & unpacking variables into structures and
            %    vice versa, according to the syntax and inputs.
            %
            %    Function features:
            %      * Pack variables to structure with enhanced field naming
            %      * Pack and update variables in existing structure
            %      * Unpack variables from structure with enhanced variable naming
            %      * Unpack only specific fields in a structure to variables
            %      * Unpack without over writing existing variables in work space
            %
            %    In addition to the obvious usage, this function could by highly useful for example in
            %    working with a function with multiple inputs. Packing variables before the call to
            %    the function, and unpacking it in the beginning of the function will make the function
            %    call shorter, more readable, and you would not have to worry about arguments order any
            %    more. Moreover you could leave the function as it is and you could pass same inputs to
            %    multiple functions, each of which will use its designated arguments placed in the
            %    structure.
            %
            %% Syntax
            %    Pack
            %      S = v2struct
            %      S = v2struct(x,y,z,...)
            %      S = v2struct(fieldNames)
            %      S = v2struct(A,B,C,..., fieldNames)
            %      S = v2struct(x,..., nameOfStruct2Update, fieldNames)
            %      v2struct
            %      v2struct(x,y,z,...)
            %      v2struct(fieldNames)
            %      v2struct(A,B,C,..., fieldNames)
            %      v2struct(x,..., nameOfStruct2Update, fieldNames)
            %
            %    Unpack
            %      v2struct(S)
            %      [a,b,c,...] = v2struct(S)
            %      v2struct(S,fieldNames)
            %      [a,b,c,...] = v2struct(S,fieldNames)
            %
            %% Inputs & Outputs
            %    Pack - inputs
            %      x,y,z,...           - any variable to pack. can be replaced by fieldNames below.
            %      nameOfStruct2Update - optional, name of structure to update if desired.
            %      fieldNames          - optional, cell array of strings, which must include a cell
            %                            with the string 'fieldNames' and must be the last input.
            %    Pack - outputs
            %      S - the packed structure. If there is no output argument then a structure named
            %          Sv2struct would be created in the caller workspace.
            %
            %    Unpack - inputs
            %      S          - name of structure to be unpacked.
            %      fieldNames - optional, cell array of strings, which must include a cell with the
            %                   string 'fieldNames' and must be the last input.
            %    Unpack - outputs
            %      a,b,c,... - variables upacked from the structure.
            %                  if there are no output arguments then variables would be created in
            %                  the caller workspace with naming according to name of inputs.
            %
            %% Examples
            %  % see 'Usage example' section below for convenient presentation of these examples.
            %
            %    % NOTE: whenever using filedNames cell array please note the following
            %    %  1. fieldNames cell array must include a cell with the string 'fieldNames'
            %    %  2. fieldNames cell array input must be the last input.
            %
            %  % Pack
            %      x = zeros(3); x2 = ones(3); y = 'Testing123'; z = cell(2,3);
            %      fieldNames1 = {'fieldNames', 'x', 'y', 'z'};
            %      fieldNames2 = {'fieldNames', 'a', 'b', 'c'};
            %      fieldNames3 = {'fieldNames', 'x'};
            %      nameOfStruct2Update = 'S';
            %
            %       % The four examples below return structure S with same values however the
            %       % structure's field names are defined differently in every syntax.
            %      % Example 1.
            %      % structure field names defined by variables names.
            %       S = v2struct(x,y,z)
            %      % Example 2.
            %      % structure field names defined according to the cell array fieldNames.
            %       % NOTE: variables with the names in fieldNames1 must exist in the caller workspace.
            %       S = v2struct(fieldNames1)
            %      % Example 3.
            %      % same as #1. but arguments are passed explicitly
            %       S = v2struct(zeros(3), 'Testing123', cell(2,3), fieldNames1)
            %      % Example 4.
            %      % field names defined by content of fieldNames2 while
            %      % the values are set according to the passed arguments. In this case the structure
            %      % S returned would be: S.a=x, S.b=y, S.c=z
            %       S = v2struct(x,y,z, fieldNames2)
            %
            %      % Example 5.
            %      % update structure S. The fields that would be updated are according to content
            %      % of fieldNames3. Note that you must pass a variable with the name
            %      % 'nameOfStruct2Update' placed before 'fieldNames3'. This variable should contain
            %      % the name of the structure you want to update as a string. Also note that if you
            %      % set an output structure name which is different than the one stated in
            %      % nameOfStruct2Update a new structure would be created and the structure that was
            %      % meant to be updated would not get updated.
            %       S.oldField = 'field to be saved for future use'
            %       S = v2struct(x2, nameOfStruct2Update, fieldNames3)
            %
            %      % Example 6.
            %      % pack all variables in caller workspace. Call without input arguments.
            %        S = v2struct
            %
            %      % The following examples return the same results as the examples above but the
            %      % structure would be returned with the default name 'Sv2struct'. Be cautious as
            %      % this might lead to overriding of arguments.
            %      % Example 7.
            %       v2struct(x,y,z)
            %      % Example 8.
            %       v2struct(fieldNames1)
            %      % Example 9.
            %       v2struct(zeros(3), 'Testing123', cell(2,3), fieldNames1)
            %      % Example 10.
            %       v2struct(x,y,z, fieldNames2)
            %      % Example 11.
            %       S.oldField = 'field to be saved for future use'
            %       v2struct(x2, nameOfStruct2Update, fieldNames3)
            %      % Example 12.
            %       v2struct
            %
            %  % Unpack
            %      clear S x x2 y z fieldNames1 fieldNames2 fieldNames3 nameOfStruct2Update
            %      S.x = zeros(3); S.y = 'Testing123'; S.z = cell(2,3);
            %      fieldNames3 = {'fieldNames','x','z'};
            %
            %      % Example 1.
            %      % This example creates or overwrites variables x, y, z in the caller with the
            %      % contents of the corresponding named fields.
            %       v2struct(S)
            %
            %      % Example 2.
            %      % This example assigns the contents of the fields of the scalar structure
            %      % S to the variables a,b,c rather than overwriting variables in the caller. If
            %      % there are fewer output variables than there are fields in S, the remaining fields
            %      % are not extracted.
            %       [a,b,c] = v2struct(S)
            %
            %      % Example 3.
            %      % This example creates or overwrites variables x and z in the caller with the
            %      % contents of the corresponding named fields.
            %       v2struct(S, fieldNames3)
            %
            %      % Example 4.
            %      % This example assigns the contents of the fields 'x' and 'z' defined by
            %      % fieldNames3 of the scalar structure S to the variables a and b rather than
            %      % overwriting variables in the caller. If there are fewer output variables than
            %      % there are fields in S, the remaining fields are not extracted.
            %       [a,b] = v2struct(S, fieldNames3)
            %
            %       % This example unpacks variables 'y' and 'z' only without overwriting variable 'x'.
            %       % NOTE the addition of the field named 'avoidOverWrite' to the structure to be
            %       % unpacked. This is mandatory in order to make this functionality work. The
            %       % contents of this field can be anything, it does not matter.
            %      S.avoidOverWrite = '';
            %      x = 'do not overwrite me';
            %      v2struct(S)
            %
            %% Usage example (includes sub-functions)
            %    1. run attached v2structDemo1.m file for on screen presentation of examples.
            %    2. run attached v2structDemo2.m file and read comments in file for a suggestion of
            %       how to use v2struct in managing input to other functions with improved usability.
            %
            %% Revision history
            %    2011-05-19, Adi N., Creation
            %    2011-05-29, Adi N., Update structure added, some documentation and demo function changes
            %    2011-06-02, Adi N., Fixed updating structure functionality
            %    2011-06-05, Adi N., Added functionality: avoid overwritring existing variables, added
            %                        unpacking examples to demo1 .m file.
            %    2011-06-30, Adi N., fieldNames usage corrected, now must include a specific string to
            %                        be triggered. Documentation enhanced. Code tweaked.
            %    2011-07-14, Adi N., Fixed bug in packing with variables only
            %    2011-08-14, Adi N., Clarified warning and error when packing/unpacking with
            %                        fieldNames.
            %    2011-09-12, Adi N., Added easy packing of all variables in caller workspace (thanks
            %                        to Vesa Lehtinen for the suggestion), fixed bug in warning
            %                        handling in packing case, edited comments.
            %
            %    Inspired by the function: mmv2truct - D.C. Hanselman, University of Maine, Orono, ME
            %    04469 4/28/99, 9/29/99, renamed 10/19/99 Mastering MATLAB 5, Prentice Hall,
            %    ISBN 0-13-858366-8
            
            
            % parse input for field names
            if isempty(varargin)
                gotCellArrayOfStrings = false;
                toUnpackRegular = false;
                toUnpackFieldNames = false;
                gotFieldNames = false;
            else
                gotCellArrayOfStrings = iscellstr(varargin{end});
                
                toUnpackRegular = (nargin == 1) && isstruct(varargin{1});
                if toUnpackRegular
                    fieldNames = fieldnames(varargin{1})';
                    nFields = length(fieldNames);
                end
                
                gotFieldNames = gotCellArrayOfStrings & any(strcmpi(varargin{end},'fieldNames'));
                if gotFieldNames
                    fieldNamesRaw = varargin{end};
                    % indices of cells with actual field names, leaving out the index to 'fieldNames' cell.
                    indFieldNames = ~strcmpi(fieldNamesRaw,'fieldNames');
                    fieldNames = fieldNamesRaw(indFieldNames);
                    nFields = length(fieldNames);
                end
                toUnpackFieldNames = (nargin == 2) && isstruct(varargin{1}) && gotFieldNames;
            end
            
            
            % Unpack
            if toUnpackRegular || toUnpackFieldNames
                
                struct = varargin{1};
                assert(isequal(length(struct),1) , 'Single input nust be a scalar structure.');
                CallerWS = evalin('caller','whos'); % arguments in caller work space
                
                % update fieldNames according to 'avoidOverWrite' flag field.
                if isfield(struct,'avoidOverWrite')
                    indFieldNames = ~ismember(fieldNames,{CallerWS(:).name,'avoidOverWrite'});
                    fieldNames = fieldNames(indFieldNames);
                    nFields = length(fieldNames);
                end
                
                if toUnpackRegular % Unpack with regular fields order
                    if nargout == 0 % assign in caller
                        for iField = 1:nFields
                            assignin('caller',fieldNames{iField},struct.(fieldNames{iField}));
                        end
                    else % dump into variables
                        for iField = 1:nargout
                            varargout{iField} = struct.(fieldNames{iField});
                        end
                    end
                    
                elseif toUnpackFieldNames % Unpack with fields according to fieldNames
                    if nargout == 0 % assign in caller, by comparing fields to fieldNames
                        for iField = 1:nFields
                            assignin('caller',fieldNames{iField},struct.(fieldNames{iField}));
                        end
                    else % dump into variables
                        assert( isequal(nFields, nargout) , ['Number of output arguments',...
                            ' does not match number of field names in cell array']);
                        for iField = 1:nFields
                            varargout{iField} = struct.(fieldNames{iField});
                        end
                    end
                end
                
                % Pack
            else
                % build cell array of input names
                CallerWS = evalin('caller','whos');
                inputNames = cell(1,nargin);
                for iArgin = 1:nargin
                    inputNames{iArgin} = inputname(iArgin);
                end
                nInputs = length(inputNames);
                
                % look for 'nameOfStruct2Update' variable and get the structure name
                if ~any(strcmpi(inputNames,'nameOfStruct2Update')) % no nameOfStruct2Update
                    nameStructArgFound = false;
                    validVarargin = varargin;
                else % nameOfStruct2Update found
                    nameStructArgFound = true;
                    nameStructArgLoc = strcmp(inputNames,'nameOfStruct2Update');
                    nameOfStruct2Update = varargin{nameStructArgLoc};
                    % valid varargin with just the inputs to pack and fieldNames if exists
                    validVarargin = varargin(~strcmpi(inputNames,'nameOfStruct2Update'));
                    % valid inputNames with just the inputs name to pack and fieldNames if exists
                    inputNames = inputNames(~strcmpi(inputNames,'nameOfStruct2Update'));
                    nInputs = length(inputNames);
                    % copy structure from caller workspace to enable its updating
                    if ismember(nameOfStruct2Update,{CallerWS(:).name}) % verify existance
                        S = evalin('caller',nameOfStruct2Update);
                    else
                        error(['Bad input. Structure named ''',nameOfStruct2Update,...
                            ''' was not found in workspace'])
                    end
                end
                
                % when there is no input or the input is only variables and perhaps
                % also nameOfStruct2Update
                if ~gotFieldNames
                    % no input, pack all of variables in caller workspace
                    if isequal(nInputs, 0)
                        for iVar = 1:length(CallerWS)
                            S.(CallerWS(iVar).name) = evalin('caller',CallerWS(iVar).name);
                        end
                        % got input, check input names and pack
                    else
                        for iInput = 1:nInputs
                            if gotCellArrayOfStrings % called with a cell array of strings
                                errMsg = sprintf(['Bad input in cell array of strings.'...
                                    '\nIf you want to pack (or unpack) using this cell array as'...
                                    ' designated names'...
                                    '\nof the structure''s fields, add a cell with the string'...
                                    ' ''fieldNames'' to it.']);
                            else
                                errMsg = sprintf(['Bad input in argument no. ', int2str(iArgin),...
                                    ' - explicit argument.\n'...
                                    'Explicit arguments can only be called along with a matching'...
                                    '\n''fieldNames'' cell array of strings.']);
                            end
                            assert( ~isempty(inputNames{iInput}), errMsg);
                            S.(inputNames{iInput}) = validVarargin{iInput};
                        end
                        
                        % issue warning for possible wrong usage when packing with an input of cell array of
                        % strings as the last input without it containing the string 'fieldNames'.
                        if gotCellArrayOfStrings
                            name = inputNames{end};
                            % input contains structure and a cell array of strings
                            if (nargin == 2) && isstruct(varargin{1})
                                msgStr = [inputNames{1},''' and ''',inputNames{2},''' were'];
                                % input contains any arguments with an implicit cell array of strings
                            else
                                msgStr = [name, ''' was'];
                            end
                            warnMsg = ['V2STRUCT - ''%s packed in the structure.'...
                                '\nTo avoid this warning do not put ''%s'' as last v2struct input.'...
                                '\nIf you want to pack (or unpack) using ''%s'' as designated names'...
                                ' of the'...
                                '\nstructure''s fields, add a cell with the string ''fieldNames'' to'...
                                ' ''%s''.'];
                            fprintf('\n')
                            warning('MATLAB:V2STRUCT:cellArrayOfStringNotFieldNames',warnMsg,msgStr,...
                                name,name,name)
                        end
                    end
                    % fieldNames cell array exists in input
                elseif gotFieldNames
                    nVarToPack = length(varargin)-1-double(nameStructArgFound);
                    if nVarToPack == 0 % no variables to pack
                        for iField = 1:nFields
                            S.(fieldNames{iField}) = evalin('caller',fieldNames{iField});
                        end
                        
                        % else - variables to pack exist
                        % check for correct number of fields vs. variables to pack
                    elseif ~isequal(nFields,nVarToPack)
                        error(['Bad input. Number of strings in fieldNames does not match',...
                            'number of input arguments for packing.'])
                    else
                        for iField = 1:nFields
                            S.(fieldNames{iField}) = validVarargin{iField};
                        end
                    end
                    
                end % if ~gotFieldNames
                
                if nargout == 0
                    assignin( 'caller', 'Sv2struct',S );
                else
                    varargout{1} = S;
                end
                
            end % if nargin
            
        end
        
        
        
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        
        % figset(hh,xlab,ylab,tt,fntsz,varargin)
        %
        % utility function which provides my favorite figure settings
        %
        % - hh: current figure handle (e.g., gca)
        % - xlab: (string) xlabel
        % - ylab: (string) ylabel
        % - tt:   (string) title
        % - fntsz: (numeric) font size (typically between 15 and 22)
        % - varargin{1}: (string) zlabel in 3D plots
        %
        % USAGE: figset(gca,xlab,ylab,tt,fntsz,varargin);
        %
        % GLC, 20 Nov 2007 -- added optional color parameter Aug 2012
        
        function figset(hh,xlab,ylab,tt,fntsz,varargin)
            col='none';
            xlabel(xlab,'fontsize',fntsz);
            ylabel(ylab,'fontsize',fntsz);
            title(tt,'fontsize',fntsz);
            set(hh,'fontsize',fntsz,'color',col);
            box off;
            if ~isempty(varargin) % 2D plot
                zlab=varargin{1};
                zlabel(zlab,'fontsize',fntsz);
            end
        end
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        % function [tout z] = gaussfilt(t,x,binwidth,causal_filter)
        %
        % filter 1D vector x with a Gaussian filter of std = binwidth (half width is
        % half of that) by convolution of x with a gaussian. It uses 'conv'.
        % Based on
        % http://www.mathworks.com/matlabcentral/newsreader/view_thread/272556
        %
        % INPUT:
        % t = time vector
        % x = data
        % binwidth = twice the halfwidth (= std of gaussian filter). NOTE: must be
        % integer > 1; width of gaussian kernel for filtering (same units as t);
        % set it to 1 to prevent filtering
        % causal_filter = flag: if ==1, causal filter used (i.e., =0 at negative
        % times).
        %
        % NOTE: a fictitious continuation of the t and x values to the left and
        % right of the data vector are added to prevent bad filtering at the
        % boundaries.
        %
        % GLC, Nov 7, 2013
        
        function  [tout z] = gaussfilt(t,x,binwidth,causal_filter)
            
            
            if size(t,1)==1 t=t'; end;
            if size(x,1)==1 x=x'; end;
            exlen=binwidth;
            x1=repmat(x(1),exlen,1);
            x2=repmat(x(end),exlen,1);
            x0=[x1 ; x ; x2];
            tbin=diff(t(1:2));
            t1=fliplr(t(1):-tbin:t(1)-(exlen-1)*tbin);
            t2=fliplr(t(end):tbin:t(end)+(exlen-1)*tbin);
            t0=[t1' ; t ; t2'];
            
            % Construct blurring window:
            windowWidth = int16(binwidth);
            halfWidth = windowWidth / 2;
            gaussFilter = gausswin(binwidth);
            if causal_filter==1; gaussFilter(1:floor(max(size(gaussFilter))/2))=0; end;
            if causal_filter==-1; gaussFilter(floor(max(size(gaussFilter))/2)+1:end)=0; end;
            gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
            
            % alpha filter -- IN PROGRESS (13 nov 2013) -- do not trust
            % bw=round(binwidth/2);
            % alphac=1/bw/bw; y=alphac^2*(1:2*bw/alphac).*exp(-alphac*(1:2*bw/alphac));
            % y=[ zeros(length(y),1)' y];
            % masky=1:floor(length(y)/binwidth):length(y);
            % y=y(masky);
            % gaussFilter=y/sum(y);
            
            % Do the blur:
            z = conv(x0, gaussFilter);
            z = z(halfWidth+exlen:end-halfWidth-exlen+1);
            % tout = t0(halfWidth+exlen:end-halfWidth-exlen+1);
            tout=exlen*tbin+t0(1:length(z));
            
            
            % plot:
            % figure(1); clf;
            % hold on;
            % plot(t,x,'k');
            % plot(t0,z, 'r-', 'linewidth', 3);
            % hold off;
        end
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        function name = getComputerName()
            % GETCOMPUTERNAME returns the name of the computer (hostname)
            % name = getComputerName()
            %
            % WARN: output string is converted to lower case
            %
            %
            % See also SYSTEM, GETENV, ISPC, ISUNIX
            %
            % m j m a r i n j (AT) y a h o o (DOT) e s
            % (c) MJMJ/2007
            % MOD: MJMJ/2013
            
            [ret, name] = system('hostname');
            
            if ret ~= 0,
                if ispc
                    name = getenv('COMPUTERNAME');
                else
                    name = getenv('HOSTNAME');
                end
            end
            ind=strfind(name,'-');
            if ~isempty(ind)
                name=name(1:ind(1)-1);
            end
            name = strtrim(lower(name));
        end
        
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        
        
        %MTIT		creates a major title in a figure with many axes
        %
        %		MTIT
        %		- creates a major title above all
        %		  axes in a figure
        %		- preserves the stack order of
        %		  the axis handles
        %
        %SYNTAX
        %-------------------------------------------------------------------------------
        %		P = MTIT(TXT,[OPT1,...,OPTn])
        %		P = MTIT(FH,TXT,[OPT1,...,OPTn])
        %
        %INPUT
        %-------------------------------------------------------------------------------
        %    FH	:	a valid figure handle		[def: gcf]
        %   TXT	:	title string
        %
        % OPT	:	argument
        % -------------------------------------------
        %  xoff	:	+/- displacement along X axis
        %  yoff	:	+/- displacement along Y axis
        %  zoff	:	+/- displacement along Z axis
        %
        %		title modifier pair(s)
        % -------------------------------------------
        %   TPx	:	TVx
        %		see: get(text) for possible
        %		     parameters/values
        %
        %OUTPUT
        %-------------------------------------------------------------------------------
        % par	:	parameter structure
        %  .pos :	position of surrounding axis
        %   .oh	:	handle of last used axis
        %   .ah :	handle of invisible surrounding axis
        %   .th :	handle of main title
        %
        %EXAMPLE
        %-------------------------------------------------------------------------------
        %	subplot(2,3,[1 3]);		title('PLOT 1');
        %	subplot(2,3,4); 		title('PLOT 2');
        %	subplot(2,3,5); 		title('PLOT 3');
        %	axes('units','inches',...
        %	     'color',[0 1 .5],...
        %	     'position',[.5 .5 2 2]);	title('PLOT 41');
        %	axes('units','inches',...
        %	     'color',[0 .5 1],...
        %	     'position',[3.5 .5 2 2]);	title('PLOT 42');
        %	shg;
        %	p=mtit('the BIG title',...
        %	     'fontsize',14,'color',[1 0 0],...
        %	     'xoff',-.1,'yoff',.025);
        % % refine title using its handle <p.th>
        %	set(p.th,'edgecolor',.5*[1 1 1]);
        
        % created:
        %	us	24-Feb-2003		/ R13
        % modified:
        %	us	24-Feb-2003		/ CSSM
        %	us	06-Apr-2003		/ TMW
        %	us	13-Nov-2009 17:38:17
        %      lev      08-Apr-2011             / Added author and license info
        
        %   Author: Urs (us) Schwarz
        %   Copyright 2009 Urs Schwarz
        %   http://www.mathworks.com/matlabcentral/fileexchange/3218-mtit-a-pedestrian-major-title-creator
        
        function	par=mtit(varargin)
            
            defunit='normalized';
            if	nargout
                par=[];
            end
            
            % check input
            if	nargin < 1
                help(mfilename);
                return;
            end
            if	isempty(get(0,'currentfigure'))
                disp('MTIT> no figure');
                return;
            end
            
            vl=true(size(varargin));
            if	ischar(varargin{1})
                vl(1)=false;
                figh=gcf;
                txt=varargin{1};
            elseif	any(ishandle(varargin{1}(:)))		&&...
                    ischar(varargin{2})
                vl(1:2)=false;
                figh=varargin{1};
                txt=varargin{2};
            else
                error('MTIT> invalid input');
            end
            vin=varargin(vl);
            [off,vout]=utils.get_off(vin{:});
            
            % find surrounding box
            ah=findall(figh,'type','axes');
            if	isempty(ah)
                disp('MTIT> no axis');
                return;
            end
            oah=ah(1);
            
            ou=get(ah,'units');
            set(ah,'units',defunit);
            ap=get(ah,'position');
            if	iscell(ap)
                ap=cell2mat(get(ah,'position'));
            end
            ap=[	min(ap(:,1)),max(ap(:,1)+ap(:,3)),...
                min(ap(:,2)),max(ap(:,2)+ap(:,4))];
            ap=[	ap(1),ap(3),...
                ap(2)-ap(1),ap(4)-ap(3)];
            
            % create axis...
            xh=axes('position',ap);
            % ...and title
            th=title(txt,vout{:});
            tp=get(th,'position');
            set(th,'position',tp+off);
            set(xh,'visible','off','hittest','on');
            set(th,'visible','on');
            
            % reset original units
            ix=find(~strcmpi(ou,defunit));
            if	~isempty(ix)
                for	i=ix(:).'
                    set(ah(i),'units',ou{i});
                end
            end
            
            % ...and axis' order
            uistack(xh,'bottom');
            axes(oah);				%#ok
            
            if	nargout
                par.pos=ap;
                par.oh=oah;
                par.ah=xh;
                par.th=th;
            end
        end
        %-------------------------------------------------------------------------------
        function	[off,vout]=get_off(varargin)
            
            % search for pairs <.off>/<value>
            
            off=zeros(1,3);
            io=0;
            for	mode={'xoff','yoff','zoff'};
                ix=strcmpi(varargin,mode);
                if	any(ix)
                    io=io+1;
                    yx=find(ix);
                    ix(yx+1)=1;
                    off(1,io)=varargin{yx(end)+1};
                    varargin=varargin(xor(ix,1));
                end
            end
            vout=varargin;
        end
        
        
        %-------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------
        
        
        function [locs pks]=peakseek(x,minpeakdist,minpeakh)
            % Alternative to the findpeaks function.  This thing runs much much faster.
            % It really leaves findpeaks in the dust.  It also can handle ties between
            % peaks.  Findpeaks just erases both in a tie.  Shame on findpeaks.
            %
            % x is a vector input (generally a timecourse)
            % minpeakdist is the minimum desired distance between peaks (optional, defaults to 1)
            % minpeakh is the minimum height of a peak (optional)
            %
            % (c) 2010
            % Peter O'Connor
            % peter<dot>ed<dot>oconnor .AT. gmail<dot>com
            
            if size(x,2)==1, x=x'; end
            
            % Find all maxima and ties
            locs=find(x(2:end-1)>=x(1:end-2) & x(2:end-1)>=x(3:end))+1;
            
            if nargin<2, minpeakdist=1; end % If no minpeakdist specified, default to 1.
            
            if nargin>2 % If there's a minpeakheight
                locs(x(locs)<=minpeakh)=[];
            end
            
            if minpeakdist>1
                while 1
                    
                    del=diff(locs)<minpeakdist;
                    
                    if ~any(del), break; end
                    
                    pks=x(locs);
                    
                    [garb mins]=min([pks(del) ; pks([false del])]); %#ok<ASGLU>
                    
                    deln=find(del);
                    
                    deln=[deln(mins==1) deln(mins==2)+1];
                    
                    locs(deln)=[];
                    
                end
            end
            
            if nargout>1,
                pks=x(locs);
            end
            
            
        end
        
        %-------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------
        
        function colors = distinguishable_colors(n_colors)
            
            if n_colors<11
                colors=[0.7000         0         0;
                    0    0.5000         0;
                    0         0    1.0000;
                    0         0    0.2069;
                    1.0000         0    0.6897;
                    0    0.3103         0;
                    1.0000    0.8276         0;
                    0    0.6207    1.0000;
                    0.6207    0.3103    0.2759;
                    0.2069    1.0000    0.7241];
            elseif n_colors>10 && n_colors<21
                colors=[    0.7000         0         0;
                    0    0.5000         0;
                    0         0    1.0000;
                    0         0    0.2069;
                    1.0000         0    0.6897;
                    0    0.3103         0;
                    1.0000    0.8276         0;
                    0    0.6207    1.0000;
                    0.6207    0.3103    0.2759;
                    0.2069    1.0000    0.7241;
                    0.4483    0.2069    0.7241;
                    0.0690    0.4828    0.5172;
                    1.0000    0.6897    1.0000;
                    0.4828    0.7241    0.1034;
                    0.8276    0.4483         0;
                    0.7586    0.6897    0.4138;
                    0.8966         0    1.0000;
                    0.1379    0.1034         0;
                    0.9310    0.0345    0.3448;
                    0.4828         0    0.3448];
                
            elseif n_colors>20 && n_colors<51
                
                
                colors=[0.7000         0         0;
                    0    0.5000         0;
                    0         0    1.0000;
                    0         0    0.2069;
                    1.0000         0    0.6897;
                    0    0.3103         0;
                    1.0000    0.8276         0;
                    0    0.6207    1.0000;
                    0.6207    0.3103    0.2759;
                    0.2069    1.0000    0.7241;
                    0.4483    0.2069    0.7241;
                    0.0690    0.4828    0.5172;
                    1.0000    0.6897    1.0000;
                    0.4828    0.7241    0.1034;
                    0.8276    0.4483         0;
                    0.7586    0.6897    0.4138;
                    0.8966         0    1.0000;
                    0.1379    0.1034         0;
                    0.9310    0.0345    0.3448;
                    0.4828         0    0.3448;
                    0.3103    0.9655    1.0000;
                    0.4828    0.4138    0.5862;
                    0.3448    0.6552    0.4483;
                    0.4138    0.3103         0;
                    0.8621    1.0000         0;
                    0.6207         0         0;
                    1.0000    0.6897    0.6897;
                    0.7931    1.0000    0.6207;
                    0    0.2759    0.6207;
                    0.9310    0.4483    1.0000;
                    0.5862    0.7931    0.9655;
                    0.9310    0.4138    0.6207;
                    0.4138    1.0000    0.4483;
                    0.5172    0.4828    0.4138;
                    1.0000    0.4828    0.3793;
                    0.1379    0.4483    1.0000;
                    0.2414         0    0.1034;
                    0.0690    0.2069    0.2759;
                    0         0    0.3793;
                    0    0.2414    0.9310;
                    0.6552    0.6207         0;
                    0.5862    0.4828    0.8621;
                    0    0.6897    0.2759;
                    0.6897    0.9310    0.7931;
                    1.0000    1.0000    0.4483;
                    0.1724    0.3103    0.2069;
                    0.6897    0.1034    0.6207;
                    1.0000    0.7586    0.3448;
                    0.4138    0.5517    0.1724;
                    0.6207         0    0.2414];
            end
            colors=colors(1:n_colors,1:3);
        end
        
        
        %-------------------------------------------------------------------------------
        %-------------------------------------------------------------------------------
        
        function c = redblue(m)
            %REDBLUE    Shades of red and blue color map
            %   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
            %   The colors begin with bright blue, range through shades of
            %   blue to white, and then through shades of red to bright red.
            %   REDBLUE, by itself, is the same length as the current figure's
            %   colormap. If no figure exists, MATLAB creates one.
            %
            %   For example, to reset the colormap of the current figure:
            %
            %             colormap(redblue)
            %
            %   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG,
            %   COLORMAP, RGBPLOT.
            
            %   Adam Auton, 9th October 2009
            
            if nargin < 1, m = size(get(gcf,'colormap'),1); end
            
            if (mod(m,2) == 0)
                % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
                m1 = m*0.5;
                r = (0:m1-1)'/max(m1-1,1);
                g = r;
                r = [r; ones(m1,1)];
                g = [g; flipud(g)];
                b = flipud(r);
            else
                % From [0 0 1] to [1 1 1] to [1 0 0];
                m1 = floor(m*0.5);
                r = (0:m1-1)'/max(m1,1);
                g = r;
                r = [r; ones(m1+1,1)];
                g = [g; 1; flipud(g)];
                b = flipud(r);
            end
            
            c = [r g b];
            
        end
        
        
        
        
    end
end