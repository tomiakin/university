%% http://www.isc.ac.uk/iscbulletin/search/catalogue/csvoutput/
% INPUT

bot_lat   = num2str(minlat);
top_lat   = num2str(maxlat);
left_lon  = num2str(minlon);
right_lon = num2str(maxlon);

min_mag=num2str(Magnitude);

%%
A=urlread(['http://isc-mirror.iris.washington.edu/cgi-bin/web-db-v4?out_format=CATQuakeML',...
    '&request=REVIEWED&searchshape=RECT',...
    '&bot_lat=',bot_lat,'&top_lat=',top_lat,'&left_lon=',left_lon,'&right_lon=',right_lon,'',...
    '&start_year=',start_year,'&start_month=',start_month,'&start_day=',start_day,'&start_time=',start_time,'',...
    '&end_year=',end_year,'&end_month=',end_month,'&end_day=',end_day,'&end_time=',end_time,'',...
    '&min_mag=',min_mag,'&req_mag_type=Any']);


%%
expr = '<event publicID=".+?>.+?</event>';
objectStrings = regexp(A,expr,'match');
Nos = length(objectStrings);
disp(['Number of Events = ',num2str(Nos)])
%%
counter = 0;
for ii = 1:Nos
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bucket = regexp(objectStrings{ii},'<description.*?>.+?</description>','match');
    try
        text1=regexp(bucket,'<text>.+?</text>','match');
        text1=text1{1}{1};
        text1=text1(length('<text>')+1:end-length('</text>'));
    catch
        text1='undefined';
    end
    
    try
        text2=regexp(bucket,'<type>.+?</type>','match');
        text2=text2{1}{1};
        text2=text2(length('<type>')+1:end-length('</type>'));
    catch
        text2='undefined';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bucket = regexp(objectStrings{ii},'<time.*?>.+?</time>','match');
    try
        text3=regexp(bucket,'<value>.+?</value>','match');
        text3=text3{1}{1};
        text3=text3(length('<value>')+1:end-length('</value>'));
        YEAR = str2double(text3(1:4));
        MONTH = str2double(text3(6:7));
        DAY = str2double(text3(9:10));
        HOUR = str2double(text3(12:13));
        MIN = str2double(text3(15:16));
        SEC = str2double(text3(18:19));
    catch
        text3='undefined';
        YEAR = [];
        MONTH = [];
        DAY = [];
        HOUR = [];
        MIN = [];
        SEC = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bucket = regexp(objectStrings{ii},'<latitude.*?>.+?</latitude>','match');
    try
        text4=regexp(bucket,'<value>.+?</value>','match');
        text4=text4{1}{1};
        text4=text4(length('<value>')+1:end-length('</value>'));
        LATITUDE = str2double(text4);
    catch
        text4='undefined';
        LATITUDE = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bucket = regexp(objectStrings{ii},'<longitude.*?>.+?</longitude>','match');
    try
        text5=regexp(bucket,'<value>.+?</value>','match');
        text5=text5{1}{1};
        text5=text5(length('<value>')+1:end-length('</value>'));
        LONGITUDE = str2double(text5);
    catch
        text5='undefined';
        LONGITUDE = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bucket = regexp(objectStrings{ii},'<depth.*?>.+?</depth>','match');
    try
        text6=regexp(bucket,'<value>.+?</value>','match');
        text6=text6{1}{1};
        text6=text6(length('<value>')+1:end-length('</value>'));
        DEPTH = str2double(text6);
    catch
        text6='undefined';
        DEPTH = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bucket = regexp(objectStrings{ii},'<mag.*?>.+?</mag>','match');
    try
        text7=regexp(bucket,'<value>.+?</value>','match');
        text7=text7{1}{1};
        text7=text7(length('<value>')+1:end-length('</value>'));
        MAGNITUDE = str2double(text7);
    catch
        text7='undefined';
        MAGNITUDE = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bucket = regexp(objectStrings{ii},'<magnitude.*?>.+?</magnitude>','match');
    number_of_magnitude_type = length(bucket);
    
    for mmm = 1 :number_of_magnitude_type 
        bucket_temp = regexp(bucket{mmm},'<mag.*?>.+?</mag>','match');
        try
            text7=regexp(bucket_temp,'<value>.+?</value>','match');
            text7=text7{1}{1};
            text7=text7(length('<value>')+1:end-length('</value>'));
            MAGNITUDE(mmm) = str2double(text7);
        catch
            text7='undefined';
            MAGNITUDE(mmm) = [];
        end
        
        try
            text8=regexp(bucket{mmm},'<type>.+?</type>','match');
            text8=text8{1};
            text8=text8(length('<type>')+1:end-length('</type>'));
            MAG_Type(mmm).M=text8;
        catch
            text8='undefined';
            MAG_Type(mmm).M=text8;
        end
        
        try
            text9=regexp(bucket{mmm},'<stationCount>.+?</stationCount>','match');
            text9=text9{1};
            text9=text9(length('<stationCount>')+1:end-length('</stationCount>'));
            STATIONS(mmm)=str2double(text9);
        catch
            text9='undefined';
            STATIONS(mmm)=NaN;
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if number_of_magnitude_type==0
    else
    counter = counter + 1;
    ii = tmp_length+counter;
    CATALOG(ii).text1 = text1;
    CATALOG(ii).text2 = text2;
    CATALOG(ii).YEAR = YEAR;
    CATALOG(ii).MONTH = MONTH;
    CATALOG(ii).DAY = DAY;
    CATALOG(ii).HOUR = HOUR;
    CATALOG(ii).MIN = MIN;
    CATALOG(ii).SEC = SEC;
    CATALOG(ii).LATITUDE = LATITUDE;
    CATALOG(ii).LONGITUDE = LONGITUDE;
    CATALOG(ii).DEPTH = DEPTH/1000;
    
    if exist('STATIONS')==0
        for mmm = 1 :number_of_magnitude_type
            MW0(mmm)=magnitude_converter(MAGNITUDE(mmm),MAG_Type(mmm).M);
        end
    else
        if isempty(STATIONS) | isnan(STATIONS)
            for mmm = 1 :number_of_magnitude_type
                MW0(mmm)=magnitude_converter(MAGNITUDE(mmm),MAG_Type(mmm).M);
            end
        else
            for mmm = 1 :number_of_magnitude_type
            MW0(mmm)=magnitude_converter(MAGNITUDE(mmm),MAG_Type(mmm).M)*STATIONS(mmm);
            ST0(mmm)=STATIONS(mmm);
            end
        end
    end
    
    if exist('STATIONS')==0
        MW0(isnan(MW0))=[];
        CATALOG(ii).MW = round(mean(MW0)*10)/10;
    else
        if isempty(STATIONS) | isnan(STATIONS)
            MW0(isnan(MW0))=[];
            CATALOG(ii).MW = round(mean(MW0)*10)/10;
        else
            MW0(isnan(MW0))=[];
            ST0(isnan(ST0))=[];
            if sum(ST0)==0
                CATALOG(ii).MW = round(mean(MW0)*10)/10;
            else
                CATALOG(ii).MW = round(sum(MW0) / sum(ST0)*10)/10;
            end
        end
    end
    clear MAGNITUDE MAG_Type STATIONS MW0 ST0
    end
end
clear A hour month Final_date ii sec Initial_date lat_e start_day              
clear Magnitude              long_e                 start_month            
clear Rectangle_of_interest  mag                    start_time             
clear day                    magType                start_year             
clear depth                  maxlat                 tmp_length             
clear end_day                maxlon                 year                   
clear end_month              minlat                 
clear end_time               minlon                 
clear end_year               minute   

clear bot_lat                   text1                     
clear DAY                       bucket                    text2                     
clear DEPTH                     bucket_temp               text3                     
clear HOUR                      counter                   text4                     
clear LATITUDE                  expr                      text5                     
clear LONGITUDE                 left_lon                  text6                     
clear MIN                       min_mag                   text7                     
clear MONTH                     mmm                       text8                     
clear Nos                       number_of_magnitude_type  text9                     
clear SEC                       objectStrings             top_lat                   
clear YEAR                      right_lon              