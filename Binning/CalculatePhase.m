function ph = CalculatePhase(peaks)
 % peaks is a binary vector (0 and 1 values only) specifying peaks
 % as 1 - phase is interpolated with 0 at peaks, up to 2*pi between
 % 22/06/2015
 % B.stemkens@umcutrecht.nl
 
 ploc=find(peaks>0.5);
 ph = peaks*0;
 if (sum(peaks<-0.5)==0)
   %%% positive peaks only
   for n=1:(length(ploc)-1)
     for m=ploc(n):(ploc(n+1)-1)
       ph(m)=(m-ploc(n))/(ploc(n+1)-ploc(n))*2*pi;
     end
   end
   % tidy up the ends
   for m=1:ploc(1)
     ph(m)=wrap((m-ploc(1))/(ploc(2)-ploc(1))*2*pi,0,2*pi);
   end
   for m=ploc(end):length(ph)
     ph(m)=wrap((m-ploc(end))/(ploc(end)-ploc(end-1))*2*pi,0,2*pi);
   end
 else
   %%% positive and negative peaks (assumes strict alternation!)
   plocm=find(peaks<-0.5);
   if (abs(length(plocm)-length(ploc))>1)
     disp('>>ERROR - positive and negative peak numbers too dissimilar');
   end
   % do beginning part (up to first positive peak)
   if (plocm(1)<ploc(1))
     % if there is a negative peak first
     for m=1:ploc(1)
       ph(m)=wrap((m-plocm(1))/(ploc(1)-plocm(1))*pi+pi,0,2*pi);
     end
     plocm=plocm(2:end);
   else
     % if there is a positive peak first
     for m=1:ploc(1)
       ph(m)=wrap((m-ploc(1))/(plocm(1)-ploc(1))*pi,0,2*pi);
     end
   end
   % do main part
   for n=1:(length(ploc)-1)
     for m=ploc(n):(plocm(n)-1)
       ph(m)=wrap((m-ploc(n))/(plocm(n)-ploc(n))*pi,0,2*pi);
     end
     for m=plocm(n):(ploc(n+1)-1)
       ph(m)=wrap((m-plocm(n))/(ploc(n+1)-plocm(n))*pi+pi,0,2*pi);
     end
   end
   % do end (after last positive peak)
   if (length(plocm)<length(ploc))
     % last peak is positive
     for m=ploc(end):length(ph)
       ph(m)=wrap((m-ploc(end))/(ploc(end)-plocm(end))*pi,0,2*pi);
     end
   else
     % last peak is negative
     for m=ploc(end):(plocm(end)-1)
       ph(m)=wrap((m-ploc(end))/(plocm(end)-ploc(end))*pi,0,2*pi);
     end
     for m=plocm(end):length(ph)
       ph(m)=wrap((m-plocm(end))/(plocm(end)-ploc(end))*pi+pi,0,2*pi);
     end
   end
 end