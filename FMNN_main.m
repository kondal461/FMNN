clc
clear all
close all

load('Iris_dataset.mat');
A = Iris_dataset;
m=size(A,1)
n=size(A,2)
V=[];
W=[];
the_ta=0.2      
y=10;           
l=1;            
k=0;            
k1=0;           
h1=0;
h2=0;

for i=1:n
    C(1,i)=A(1,i);
    V=[V A(1,i)];
    W=[W A(1,i)];
end
BV=[V];
BW=[W];
BV(1,n+1)=l;            
BW(1,n+1)=l;            
C(1,n+1)=l;             
C(1,n+1)=l;              
for h=2:m    
   for i=1:n
        C(h,i)=A(h,i);
   end
    kt_tm=0;  
    max_bj=0; 
    max_bj1=0;
    for j=1:l
        bj=0;
        for i=1:n 
            aw=(A(h,i)-BW(j,i))*y;
            av=(BV(j,i)-A(h,i))*y;
            if aw > 1
                aw=1;
            else
                if aw < 0
                    aw = 0;
                end
            end
            if av > 1
                av=1;
            else
               if av < 0
                  av = 0;
               end
            end
            bj=bj+(1-aw-av);
        end
        bj=bj/(n);
       
        th=0;th1=0;th2=0;
        max_th=0;
        for i1=1:n 
            th1=BW(j,i1);
            if th1 < A(h,i1) 
               th1=A(h,i1);  
            end
            th2=BV(j,i1);
            if th2>A(h,i1) 
               th2=A(h,i1);  
            end
            th=th+(th1-th2);
        end
        th=th/(n);
        if  th<=the_ta
          if max_bj < bj
             max_bj=bj;
             k=j;
             kt_tm=1;
          end                    
        end               
    end
    max_bj;
    if (kt_tm==1)
        C(h,n+1)=k;
         for i2=1:n
        	if BV(k,i2)>A(h,i2) 
            	BV(k,i2)=A(h,i2); 
            end
            if BW(k,i2)<A(h,i2) 
            	BW(k,i2)=A(h,i2); 
            end          
       end 
        for j=1:l 
            if (j ~= k)&(BV(j,n-1)~=BV(k,n-1))               
                anpha_o=1;
                anpha_n=1;
                delta=0;
                xetmin=0;
                for i3=1:n 
                    if (delta>=0)  
                        if (BV(j,i3)< BV(k,i3))&(BV(k,i3)<BW(j,i3))&(BW(j,i3)< BW(k,i3))
                           if (anpha_o > BW(j,i3)- BV(k,i3))
                                anpha_n = BW(j,i3)- BV(k,i3);
                           end
                           if (anpha_o < BW(j,i3)- BV(k,i3))
                            anpha_n = anpha_o;
                           end
                        end
                        if (BV(k,i3) < BV(j,i3))&(BV(j,i3)< BW(k,i3))&(BW(k,i3)< BW(j,i3))
                            if (anpha_o > BW(k,i3)- BV(j,i3))
                                anpha_n = BW(k,i3)- BV(j,i3);
                            end
                            if (anpha_o < BW(k,i3)- BV(j,i3))
                                anpha_n = anpha_o;
                            end 
                        end
                        if (BV(j,i3) == BV(k,i3))&(BV(k,i3)<BW(j,i3))&(BW(j,i3)< BW(k,i3))
                           if (BW(j,i3)-BV(k,i3)<BW(k,i3)-BV(j,i3))
                               xetmin = BW(j,i3)-BV(k,i3);
                           else
                               xetmin = BW(k,i3)-BV(j,i3)
                           end
                           if (anpha_o > xetmin)
                                anpha_n = xetmin;
                           end
                           if (anpha_o < xetmin)
                            anpha_n = anpha_o;
                           end
                        end                        
                        if (BV(j,i3) < BV(k,i3))&(BV(k,i3)<BW(j,i3))&(BW(j,i3)== BW(k,i3))
                           if (BW(j,i3)-BV(k,i3)<BW(k,i3)-BV(j,i3))
                               xetmin = BW(j,i3)-BV(k,i3);
                           else
                               xetmin = BW(k,i3)-BV(j,i3)
                           end
                           if (anpha_o > xetmin)
                                anpha_n = xetmin;
                           end
                           if (anpha_o < xetmin)
                            anpha_n = anpha_o;
                           end
                        end                         
                        if (BV(k,i3) == BV(j,i3))&(BV(j,i3)< BW(k,i3))&(BW(k,i3)< BW(j,i3))
                           if (BW(j,i3)-BV(k,i3)<BW(k,i3)-BV(j,i3))
                               xetmin = BW(j,i3)-BV(k,i3);
                           else
                               xetmin = BW(k,i3)-BV(j,i3)
                           end
                           if (anpha_o > xetmin)
                                anpha_n = xetmin;
                           end
                           if (anpha_o < xetmin)
                            anpha_n = anpha_o;
                           end 
                        end                        
                        if (BV(k,i3) < BV(j,i3))&(BV(j,i3)< BW(k,i3))&(BW(k,i3)== BW(j,i3))
                           if (BW(j,i3)-BV(k,i3)<BW(k,i3)-BV(j,i3))
                               xetmin = BW(j,i3)-BV(k,i3);
                           else
                               xetmin = BW(k,i3)-BV(j,i3)
                           end
                           if (anpha_o > xetmin)
                                anpha_n = xetmin;
                           end
                           if (anpha_o < xetmin)
                            anpha_n = anpha_o;
                           end 
                        end                         
                        if (BV(j,i3)<BV(k,i3))&(BV(k,i3)<=BW(k,i3))&(BW(k,i3)<BW(j,i3))
                            min = BW(k,i3)-BV(j,i3);
                            if (min > (BW(j,i3)-BV(k,i3)))
                                min = BW(j,i3)-BV(k,i3);
                            end                         
                            if (anpha_o > min)
                            	anpha_n = min;
                            end
                            if (anpha_o < min)
                                anpha_n = anpha_o;
                            end
                        end 
                        if (BV(k,i3)<BV(j,i3))&(BV(j,i3)<=BW(j,i3))&(BW(j,i3)<BW(k,i3))
                            min = BW(j,i3)-BV(k,i3)                        
                            if (min>(BW(k,i3)-BV(j,i3)))
                                min = BW(k,i3)-BV(j,i3)
                            end
                            if (anpha_o > min)
                                anpha_n = min;                                
                            end
                            if (anpha_o < min)
                                anpha_n = anpha_o;
                            end      
                        end
                        if (BV(k,i3) == BV(j,i3))&(BV(j,i3)< BW(k,i3))&(BW(k,i3)== BW(j,i3))
                           if (anpha_o > BW(k,i3)-BV(j,i3))
                                anpha_n = BW(k,i3)-BV(j,i3);
                           end
                           if (anpha_o < BW(k,i3)-BV(j,i3))
                            anpha_n = anpha_o;
                           end 
                        end                                                 
                        anpha_n;
                        anpha_o   ;                    
                        if ((anpha_o - anpha_n)>0)
                            delta=i3;
                            anpha_n=1;
                            anpha_o=anpha_n;                 
                        else
                            delta=-1;
                        end 
                    end  
             end
             if (delta>0)               
                if ((BV(j,delta)<BV(k,delta)) & (BV(k,delta)<BW(j,delta)) & (BW(j,delta)<BW(k,delta))) 
                         BW(j,delta)= (BW(j,delta)+BV(k,delta))/2                    
                         BV(k,delta)= BW(j,delta)
                end
                if ((BV(k,delta)<BV(j,delta))&(BV(j,delta)<BW(k,delta))&(BW(k,delta)<BW(j,delta)))  
                        BW(k,delta)= (BW(k,delta)+BV(j,delta))/2                    
                        BV(j,delta)= BW(k,delta)
                end
                if ((BV(j,delta) == BV(k,delta)) & (BV(k,delta)<BW(j,delta)) & (BW(j,delta)<BW(k,delta)))  
                         BW(j,delta)= (BW(j,delta)+BV(k,delta))/2                    
                         BV(k,delta)= BW(j,delta)
                end                
                if ((BV(j,delta) < BV(k,delta)) & (BV(k,delta)<BW(j,delta)) & (BW(j,delta)==BW(k,delta))) 
                         BW(j,delta)= (BW(j,delta)+BV(k,delta))/2                    
                         BV(k,delta)= BW(j,delta)
                end                
                if ((BV(k,delta)==BV(j,delta))&(BV(j,delta)<BW(k,delta))&(BW(k,delta)<BW(j,delta)))   
                        BW(k,delta)= (BW(k,delta)+BV(j,delta))/2                    
                        BV(j,delta)= BW(k,delta)
                end                
                if ((BV(k,delta)<BV(j,delta))&(BV(j,delta)<BW(k,delta))&(BW(k,delta)==BW(j,delta))) 
                        BW(k,delta)= (BW(k,delta)+BV(j,delta))/2                    
                        BV(j,delta)= BW(k,delta)
                end                
                if ((BV(j,delta)<BV(k,delta))&(BV(k,delta)<=BW(k,delta))&(BW(k,delta)<BW(j,delta)))
                   if (BW(k,delta)-BV(j,delta)<BW(j,delta)-BV(k,delta))
                            BV(j,delta)=BW(k,delta)
                    end
                    if (BW(k,delta)-BV(j,delta)>BW(j,delta)-BV(k,delta))
                            BW(j,delta)=BV(j,delta)
                    end
                end
                if ((BV(k,delta)<BV(j,delta))&(BV(j,delta)<=BW(j,delta))&(BW(j,delta)<BW(k,delta)))
                    if (BW(k,delta)-BV(j,delta)<BW(j,delta)-BV(k,delta))
                            BW(k,delta)=BV(j,delta)
                    end
                    if (BW(k,delta)-BV(j,delta)>BW(j,delta)-BV(k,delta))
                            BV(k,delta)=BW(j,delta)
                    end   
                end 
                if(BV(k,i3) == BV(j,i3))&(BV(j,i3)< BW(k,i3))&(BW(k,i3)== BW(j,i3))
                    BV(k,delta)=(BW(j,delta)+BV(k,delta))/2                    
                    BW(j,delta)=BV(k,delta)
                end
                if(BV(j,i3) == BV(k,i3))&(BV(k,i3)< BW(j,i3))&(BW(j,i3)== BW(k,i3))
                    BV(j,delta)=(BW(k,delta)+BV(j,delta))/2                    
                    BW(k,delta)=BV(j,delta)
                end                
            end    
           end 
       end        
    else
        l=l+1;
        V=[];
        W=[];
        for i_new=1:n
            V=[V A(h,i_new)];
            W=[W A(h,i_new)];     
        end
        C(h,n+1)=l;
        V=[V 0];
        W=[W 0];
        BV=[BV;V];
        BW=[BW;W];
        BV(l,n+1)=l;
        BW(l,n+1)=l;        
    end
end
for j=1:l    
    for i=1:n
        BV(j,i)=BV(j,i)*1.05;
        BW(j,i)=BW(j,i)*0.95;
    end
end
for h=1:m    
    max_bj=0; 
    for j=1:l
        bj=0;
        for i=1:n 
            aw=(A(h,i)-BW(j,i))*y;
            av=(BV(j,i)-A(h,i))*y;
            if aw > 1
                aw=1;
            else
                if aw < 0
                    aw = 0;
                end
            end
            if av > 1
                av=1;
            else
               if av < 0
                  av = 0;
               end
            end
            bj=bj+(1-aw-av);
        end
        bj=bj/(n);
        if max_bj < bj
             max_bj=bj;
             k=j;
        end                    
    end
    if  max_bj==1
        h1=h1+1;
        for i=1:n
            C1(h1,i)=A(h,i);
        end
        C1(h1,n+1)=k;
    else
        h2=h2+1;
        for i=1:n
            C2(h2,i)=A(h,i);
        end 
        C2(h2,i)=0;
    end  
end 