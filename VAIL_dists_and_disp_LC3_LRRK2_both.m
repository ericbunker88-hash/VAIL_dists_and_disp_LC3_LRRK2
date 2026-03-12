clear,clc
curdir=pwd;
addpath('/Users/bunkeren/Documents/MATLAB/bfmatlab 2')
addpath /Users/bunkeren/Documents/MATLAB/functions
addpath('/data/bunkeren//bfmatlab 2')
addpath /data/bunkeren/functions/
warning off


wavelengths=2;
zstacks=1;
% timepoints=1;

wb=waitbar(0);


alldir=dir('2025*');
dirdir=alldir(cell2mat({alldir.isdir}));
dirnames={dirdir(1:end).name};
dirnames2=dirnames(1:end);


for y=5%2:7
    Row=char(y+64);
    for x=4%1:12
        Column=int2strz(x,2);
        


        ti=1;
        % f=1;
        for f=1%:numel(dirnames2)
            folder=dirnames2{f};
            nddir=dir([folder,'/Well',Row,Column,'*']);
            % nddir=dir(['Well',Row,Column,'*']);
        try
            % reader=bfGetReader([nddir(1).name]);
            reader=bfGetReader([folder,'/',nddir(1).name]);
        catch
            continue
        end
        
        sites=reader.getSeriesCount;
        timepoints=reader.getImageCount/wavelengths/zstacks;
        xsize=reader.getSizeX;
        ysize=reader.getSizeY;
        
        for t=timepoints
        for s=5%1:sites
%             reader.setSeries(s-1)
%             imdir=dir([Row,Column,'_s',int2str(s),'*']);
%             if numel(imdir)==0,continue,end
%             ti=1;
            
                
                tp=t+f-1;

                wb=waitbar(t/timepoints,wb,[Row,Column,' s',int2str(s)]);
                curcell=1;clear ratc rfpc bfpc
% figure                
                for z=1:zstacks
                    for c=1:wavelengths

                        if wavelengths==1
                            reader.setSeries(s-1)
                            I1(:,:,c)=bfGetPlane(reader,t,1,1,xsize,ysize);
                        else
                            imindex1=[c+(z-1)*wavelengths+(s-1)*zstacks*wavelengths+(t-1)*sites*zstacks*wavelengths];
                            s1=floor(imindex1/(reader.getImageCount+0.01));
                            reader.setSeries(s1)
                            imindex=imindex1-((reader.getImageCount)*s1);
                            I1(:,:,c)=bfGetPlane(reader,imindex,1,1,xsize,ysize);
                        end
                         
                        imsort=sort(reshape(I1(:,:,c),numel(I1(:,:,c)),1));
                        back=mode(imsort(1:numel(imsort)*0.5));
                        I2(:,:,c,z)=I1(:,:,c)-back;
% subplot(2,3,(c-1)*3+z),imshow(imadjust(I2(:,:,c,z)))
                    end
                end
                I3=max(I2,[],4);

                % dap=I3(:,:,3);%
                gfp=I3(:,:,1);%mCh-LRRK2 (dox-inducible)
                rfp=I3(:,:,2);%is actually GFP-LC3B
%                 ubc=I3(:,:,1);%640-FK2

% 
%                 bwdr=im2bw(dap,2000/65535);
%                 dapb=imgaussfilt(dap,50);
%                 dape=uint16(double(dap)./double(dapb)*1000);
%                 bwde=im2bw(dape,1250/65535);
%                 bwd1=bwdr.*bwde;%
%                 bwd2=bwareaopen(bwd1,5);
%                 bwd3=imclose(bwd2,strel('disk',3));
%                 bwd4=imfill(bwd3,'holes');
%                 bwd5=bwareaopen(bwd4,1500);
%                 bwnuc=bwd5;
                
                
                bwr1=im2bw(rfp,200/65535);
%                 bwr2=imclose(bwr1,strel('disk',5));
                bwr3=bwareaopen(bwr1,5e3);                
                bwr=bwr3;

                if 0
                    bwg1=im2bw(gfp,50/65535);
                    bwg2=bwareaopen(bwg1,5);
                    bwg3=imclose(bwg2,strel('disk',50));
                    
                    bwg4=imfill(bwg3,'holes');
                    bwg4b=bwareaopen(bwg4,100);
                    bwg5=imdilate(bwg4b,strel('disk',50));
                    bwg=bwg5;
                    bwcells0=bwr.*bwg;
                    bwcells=logical(bwareaopen(bwcells0,5e3));
                else
                    bwcells=logical(bwr);
                end
                if nnz(bwcells)<1000,continue,end

                rfpb2=imgaussfilt(rfp,25);
                rfpe=uint16(double(rfp)./double(rfpb2)*1000);

                bwre=im2bw(rfpe,2500/65535);
                bwrp1=bwre.*im2bw(rfp,1e4/65535).*bwcells;
                bwrp2=bwareaopen(bwrp1,5);

                ccrp=bwconncomp(bwrp2);
                rp_circ=regionprops(ccrp,'Circularity');

                circthresh=.25;
                longarray=find([rp_circ.Circularity]<circthresh);
                bwrp3=bwrp2;
                for lng=longarray
                    bwrp3(ccrp.PixelIdxList{lng})=0;
                end


                % bwrp3=bwrp2;
                bwrpunct=logical(bwrp3);
                ccpunct=bwconncomp(bwrpunct);
                numpunct(y,x,t,s)=ccpunct.NumObjects;
                punctsize{y,x,t,s}=[regionprops(ccpunct,'Area').Area];
                
%                 circ{y,x,s}=rp_circ.Circularity;


 if 1
                
                %single cell segmentation
                % imblur=imgaussfilt(rfp/2+dap,20);
                imblur=imgaussfilt(rfp,25);
                mounts=-double(imblur);
%                 mounts(bwnuc)=-Inf;
                ws=watershed(mounts);
                bwcells2=bwcells;
                bwcells2(ws==0)=0;
                bwcells3=bwareaopen(bwcells2,1e4);
                cc_cells=bwconncomp(bwcells3);
                
                if nnz(bwcells3)<1000,continue,end

                bwpos=zeros(size(bwcells));
                posnum=1;

                rdub=double(rfp);%udub=double(ubc);
                clear distcells tempdist udistcells
                gtemp2=zeros(size(bwpos));
                for c=1:cc_cells.NumObjects
                    
                    % if nnz(bwrpunct(cc_cells.PixelIdxList{c}))<150,continue,end
                    temp=zeros(size(bwpos));
                    temp(cc_cells.PixelIdxList{c})=1;
                    temp=logical(imfill(temp,'holes'));
                    if 1,if sum(sum(gfp(temp)))/nnz(temp)<12,continue,end,end
                    bwpos(temp)=1;
%                     gtemp(c)=sum(sum(gfp(temp)))/nnz(temp);
%                     if gtemp(c)>30,gtemp2(temp)=1;end
                    rtemp=rdub.*logical(temp)-median(rdub(logical(temp)));
                    rtemp(rtemp<0)=0;
                    temp=logical(temp);

                    tempR=corrcoef(double(gfp(temp)),double(rfp(temp)));
                    r_gr(posnum)=tempR(2);
%                     tempR=corrcoef(double(ubc(temp)),double(rfp(temp)));
%                     r_ur(posnum)=tempR(2);
%                     tempR=corrcoef(double(gfp(temp)),double(ubc(temp)));
%                     r_gu(posnum)=tempR(2);

                    props = regionprops(temp, rtemp, 'WeightedCentroid','MajorAxisLength');
% 
% 
                    rdubt=rdub(cc_cells.PixelIdxList{c});
%                     udubt=udub(cc_cells.PixelIdxList{c});
                    ypx=cc_cells.PixelIdxList{c}-floor(cc_cells.PixelIdxList{c}/size(rdub,1))*size(rdub,1);
                    xpx=floor(cc_cells.PixelIdxList{c}/size(rdub,1));
                    xdist=xpx-props.WeightedCentroid(1);
                    ydist=ypx-props.WeightedCentroid(2);
                    eud=sqrt(xdist.^2+ydist.^2);
                    tempdisp(c)=mean([eud./props.MajorAxisLength].*...
                        [rdub(cc_cells.PixelIdxList{c})./mean(rdub(cc_cells.PixelIdxList{c}))]);
                    
                    for i=1:76
                        distcells(i,posnum)=nanmean(rdubt(ceil(eud)==i));
%                         udistcells(i,posnum)=nanmean(udubt(ceil(eud)==i));
                    end
                    
                    
                    posnum=posnum+1;
                end
                blah
                if posnum==1,continue,end
                R_gr(y,x,t,s,f)=nanmedian(r_gr);
%                 R_ur(y,x,s)=nanmedian(r_ur);
%                 R_gu(y,x,s)=nanmedian(r_gu);
                disp{y,x,t,s,f}=tempdisp;
                distcell{y,x,t,s,f}=distcells;
%                 udistcell{y,x,s}=udistcells;
 end
                
        end
            end
%         end
        clear reader
%         save([Row,Column])
%         clear p4pix p7pix
        end
    end
end


clear reader
close(wb)
save(['dists_',int2str(clock)])

%%

f=1;
clear array array2 uarray2 uarray

for y=2:7
    for x=2:8
        for t=1:timepoints
            
        consol=[];consol_u=[];consolrat=[];

        for s=1:sites
            clear temp* rat

            temp=distcell{y,x,t,s,f};
%             utemp=udistcell{y,x,s};
            if numel(temp)==0,continue,end
            for i=1:size(temp,2)
                temp2(:,i)=temp(:,i)/mean(temp(:,i));
%                 utemp2(:,i)=utemp(:,i)/mean(utemp(:,i));
                % rat(:,i)=temp(2:33,i)./shiftdim([1:32]*.16,1);
                rat(i)=nanmean(temp(1:7,i))./nanmean(temp(26:32,i));
            end
            

            array(:,s)=nanmedian(temp2,2);
            consol=[consol temp2];
            consolrat=[consolrat rat];
%             uarray(:,s)=nanmedian(utemp2,2);
%             consol_u=[consol utemp2];
        end
%         if numel(temp)==0,continue,end
%         array2(:,y,x)=nanmedian(array,2);
        array2(:,y,x,t)=nanmedian(consol,2);
        % array2(:,y,x,t)=nanmedian(consolrat,2);

        % maxarray(y,x,t)=max(nanmedian(consol,2));
        meanrat(y,x,t)=nanmean(consolrat);
%         uarray2(:,y,x)=nanmedian(uarray,2);
%         uarray2(:,y,x)=nanmedian(consol_u,2);
        end
    end
end

array2=meanrat;
% array2=nanmedian(R_gr(:,:,:,:,f),4);
xs=[2:8];
of=0.145;
clear ave err xax
for y=1:2
    for x=xs
        for t=1:timepoints
            ave(t,x-xs(1)+1,y)=nanmean(array2((y-1)*3+2:y*3+1,x,t),1);
            err(t,x-xs(1)+1,y)=nanstd(array2((y-1)*3+2:y*3+1,x,t),0,1);
            
            xax(t,x-xs(1)+1,y)=(t-1)*.5+.25;

        end
    end
end


disp_array=[ave(:,:,1),ave(:,:,2)];

yarray={'WT','Rab3A','Rab3B','Rab3C','Rab3D','Rab8A','Rab8B',...
    'WT+MLi-2','Rab10','Rab12','Rab29','Rab35','Rab38','Rab43'};

my_cm=jet(100);

figure
imagesc(shiftdim(disp_array(:,:),1),[.5 2.5])
set(gcf,'colormap',my_cm)
set(gca,'ytick',[1:14],'yticklabel',yarray,...
    'xtick',1:2:13,'xticklabel',[0:6],...
    'tickdir','out')
cb=colorbar;

set(findall(gcf,'-property','fontsize'),'fontsize',20)
title('Monensin')

%%
load('dists_2025     4     8     3    56     8.mat')
f=1;
clear array array2 uarray2 uarray

for y=2:7
    for x=2:8
        for t=1:timepoints
            
        consol=[];consol_u=[];consolrat=[];

        for s=1:sites
            clear temp* rat

            temp=distcell{y,x,t,s,f};
%             utemp=udistcell{y,x,s};
            if numel(temp)==0,continue,end
            for i=1:size(temp,2)
                temp2(:,i)=temp(:,i)/mean(temp(:,i));
%                 utemp2(:,i)=utemp(:,i)/mean(utemp(:,i));
                % rat(:,i)=temp(2:33,i)./shiftdim([1:32]*.16,1);
                rat(i)=nanmean(temp(1:7,i))./nanmean(temp(26:32,i));
            end
            

            array(:,s)=nanmedian(temp2,2);
            consol=[consol temp2];
            consolrat=[consolrat rat];
%             uarray(:,s)=nanmedian(utemp2,2);
%             consol_u=[consol utemp2];
        end
%         if numel(temp)==0,continue,end
%         array2(:,y,x)=nanmedian(array,2);
        array2(:,y,x,t)=nanmedian(consol,2);
        % array2(:,y,x,t)=nanmedian(consolrat,2);

        % maxarray(y,x,t)=max(nanmedian(consol,2));
        meanrat(y,x,t)=nanmean(consolrat);
%         uarray2(:,y,x)=nanmedian(uarray,2);
%         uarray2(:,y,x)=nanmedian(consol_u,2);
        end
    end
end

array2=meanrat;
% array2=nanmedian(R_gr(:,:,:,:,f),4);
xs=[2:8];
of=0.145;
clear ave err xax
for y=1:2
    for x=xs
        for t=1:timepoints
            ave(t,x-xs(1)+1,y)=nanmean(array2((y-1)*3+2:y*3+1,x,t),1);
            err(t,x-xs(1)+1,y)=nanstd(array2((y-1)*3+2:y*3+1,x,t),0,1);
            
            xax(t,x-xs(1)+1,y)=(t-1)*.5+.25;

        end
    end
end


disp_array=[ave(:,:,1),ave(:,:,2)];
reorg_array=[1,8,2:6,9:13];
disp_array2=disp_array(:,reorg_array);

yarray={'WT','Rab3A','Rab3B','Rab3C','Rab3D','Rab8A','Rab8B',...
    'WT+MLi-2','Rab10','Rab12','Rab29','Rab35','Rab38','Rab43'};

my_cm=jet(100);

figure
imagesc(shiftdim(disp_array2(:,:),1),[.5 2.5])
set(gcf,'colormap',my_cm)
set(gca,'ytick',[1:14],'yticklabel',yarray(reorg_array),...
    'xtick',1:2:13,'xticklabel',[0:6],...
    'tickdir','out')
cb=colorbar;
xlabel('Time (hr)')
set(findall(gcf,'-property','fontsize'),'fontsize',20)
title('Monensin')
%%
load('dists_2025     4     7    18    42    23.mat')
f=2;
clear array array2 uarray2 uarray

for y=2:7
    for x=2:8
        for t=1:timepoints
            
        consol=[];consol_u=[];consolrat=[];

        for s=1:sites
            clear temp* rat

            temp=distcell{y,x,t,s,f};
%             utemp=udistcell{y,x,s};
            if numel(temp)==0,continue,end
            for i=1:size(temp,2)
                temp2(:,i)=temp(:,i)/mean(temp(:,i));
%                 utemp2(:,i)=utemp(:,i)/mean(utemp(:,i));
                % rat(:,i)=temp(2:33,i)./shiftdim([1:32]*.16,1);
                rat(i)=nanmean(temp(1:7,i))./nanmean(temp(26:32,i));
            end
            

            array(:,s)=nanmedian(temp2,2);
            consol=[consol temp2];
            consolrat=[consolrat rat];
%             uarray(:,s)=nanmedian(utemp2,2);
%             consol_u=[consol utemp2];
        end
%         if numel(temp)==0,continue,end
%         array2(:,y,x)=nanmedian(array,2);
        array2(:,y,x,t)=nanmedian(consol,2);
        % array2(:,y,x,t)=nanmedian(consolrat,2);

        % maxarray(y,x,t)=max(nanmedian(consol,2));
        meanrat(y,x,t)=nanmean(consolrat);
%         uarray2(:,y,x)=nanmedian(uarray,2);
%         uarray2(:,y,x)=nanmedian(consol_u,2);
        end
    end
end

array2=meanrat;
% array2=nanmedian(R_gr(:,:,:,:,f),4);
xs=[2:8];
of=0.145;
clear ave err xax
for y=1:2
    for x=xs
        for t=1:timepoints
            ave(t,x-xs(1)+1,y)=nanmean(array2((y-1)*3+2:y*3+1,x,t),1);
            err(t,x-xs(1)+1,y)=nanstd(array2((y-1)*3+2:y*3+1,x,t),0,1);
            
            xax(t,x-xs(1)+1,y)=(t-1)*.5+.25;

        end
    end
end


disp_array=[ave(:,:,1),ave(:,:,2)];
reorg_array=[1,8,2:6,9:14];
disp_array2=disp_array(:,reorg_array);

yarray={'WT','Rab3A','Rab3B','Rab3C','Rab3D','Rab8A','Rab8B',...
    'WT+MLi-2','Rab10','Rab12','Rab29','Rab35','Rab38','Rab43'};

my_cm=jet(100);

figure
imagesc(shiftdim(disp_array2(:,:),1),[.5 4])
set(gcf,'colormap',my_cm)
set(gca,'ytick',[1:14],'yticklabel',yarray(reorg_array),...
    'xtick',1:2:13,'xticklabel',[0:6],...
    'tickdir','out')
cb=colorbar;
set(findall(gcf,'-property','fontsize'),'fontsize',20)
title('diABZI')

%%
array2=nanmean(R_gr(:,:,:,:,f),4);

xs=[2:8];
of=0.145;
clear ave err xax
for y=1:2
    for x=xs
        for t=1:timepoints
            ave(t,x-xs(1)+1,y)=nanmean(array2((y-1)*3+2:y*3+1,x,t),1);
            err(t,x-xs(1)+1,y)=nanstd(array2((y-1)*3+2:y*3+1,x,t),0,1);
            
            xax(t,x-xs(1)+1,y)=(t-1)*.5+.25;

        end
    end
end


disp_array=[ave(:,:,1),ave(:,:,2)];



my_cm=jet(100);

figure
imagesc(shiftdim(disp_array(:,:),1),[.25 .75])
set(gcf,'colormap',my_cm)

%%

figure
errorbar(xax(:,:,1),ave(:,:,1),err(:,:,1))
hold all
errorbar(xax(:,:,1),ave(:,:,2),err(:,:,2))












