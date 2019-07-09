function [fit_poly,ADC01d,KUR1d,KUR01d,S01d,Cfit1d,errfit,ndeg,ADC0_err1d,KUR_err1d] ...
    = process_signal_POLY(data1d,bvalue,bmin,bmax)
  
  npixel = size(data1d,1);
  
  nb = length(bvalue);

  C0 =zeros(npixel,1);
  C2 =zeros(npixel,1); 
  C4 =zeros(npixel,1);
 

  ADC01d =zeros(npixel,1);
  KUR1d =zeros(npixel,1);

  fit_poly = zeros(npixel,length(bvalue));
  errfit = zeros(npixel);

  [ind_b] = find(bvalue>=bmin & bvalue <= bmax); 
  
  for ip = 1:npixel
      
    found = 0;
    Cfit_old = polyfit(bvalue(ind_b),log(data1d(ip,ind_b)),1);
    ADC0_old = -Cfit_old(1);
    KUR_old = 0;
    ndeg = -1;
    
    fit_poly = abs(polyval(Cfit_old,bvalue(ind_b)) - log(data1d(ip,ind_b)));
    
    if (max(fit_poly)<=1e-3*max(abs(log(data1d(ip,ind_b)))))
        found = 1;
        ADC0_use = ADC0_old;
        KUR_use = KUR_old;
        Cfit_use = Cfit_old;
        ndeg = 1;
        ADC0_err_use = 0;
        KUR_err_use = 0;
    end
          
    for pn = 2:length(ind_b)-1
      
      if (found == 0)     
      
        Cfit = polyfit(bvalue(ind_b),log(data1d(ip,ind_b)),pn);
        ADC0 = -Cfit(pn);
        
        KUR = Cfit(pn-1)/ADC0^2*6;
        
        ADC0_err = abs(ADC0-ADC0_old)./ADC0_old;
     
        
        KUR_err = min(abs(KUR-KUR_old),abs(KUR-KUR_old)./abs(KUR_old));
        
    
        if ((abs(ADC0-ADC0_old) <= 1e-6 | abs(ADC0-ADC0_old) <= max(0.05*ADC0_old)) & ...
              (KUR_old <= 0.15 | ((abs(KUR-KUR_old) < abs(0.15*KUR_old)) | abs(KUR-KUR_old) < 0.15)))
          found = 1;
          ADC0_use = ADC0_old;
          KUR_use = KUR_old;
          Cfit_use = Cfit_old;
          ndeg = pn-1;
          ADC0_err_use = ADC0_err;
          KUR_err_use = KUR_err;
        else
          ADC0_old = ADC0;
          KUR_old = KUR;
          Cfit_old = Cfit;
        end  
      end
    
    end
    
    if (found == 0)
      
      disp(['warning, kur may not be accurate']);
      if (exist('ADC0'))
        abs(ADC0-ADC0_old) <= max(0.01*ADC0_old) 
        (KUR_old <= 0.05 | ((abs(KUR-KUR_old) < abs(0.05*KUR_old)) | abs(KUR-KUR_old) < 0.05))
                
        ADC0
        ADC0_old
        KUR
        KUR_old
        ADC0_err
        KUR_err
        Cfit
        KUR1d(ip,1) = nan;
        KUR01d(ip,1) = nan;
        ADC01d(ip,1) = ADC0_old;
        S01d(ip,1) = exp(Cfit(end));
        Cfit1d{ip} = Cfit;
        fit_poly(ip,:) = exp(polyval(Cfit1d{ip},bvalue));
        errfit(ip,1) = sqrt(sum((data1d(ip,ind_b)-fit_poly(ip,ind_b)).^2))/length(ind_b);
        ADC0_err1d(ip,1) = nan;
        KUR_err1d(ip,1) = nan;
        %stop
      else
        KUR1d(ip,1) = nan;
        KUR01d(ip,1) = nan;
        ADC01d(ip,1) = nan;
        S01d(ip,1) = nan;
        Cfit1d{ip} = nan;
        fit_poly(ip,:) = nan;
        errfit(ip,1) = nan;
        ADC0_err1d(ip,1) = nan;
        KUR_err1d(ip,1) = nan;
      end
    else
      KUR1d(ip,1) = KUR_use;
      KUR01d(ip,1) = KUR_use;
      ADC01d(ip,1) = ADC0_use;
      S01d(ip,1) = exp(Cfit_use(end));  
      Cfit1d{ip} = Cfit_use;  
      fit_poly(ip,:) = exp(polyval(Cfit1d{ip},bvalue));
      errfit(ip,1) = sqrt(sum((data1d(ip,ind_b)-fit_poly(ip,ind_b)).^2))/length(ind_b);
      ADC0_err1d(ip,1) = ADC0_err_use;
      KUR_err1d(ip,1) = KUR_err_use;
    end
    
    
  end
  

      