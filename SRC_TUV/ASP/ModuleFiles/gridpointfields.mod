  �2  c   k820309              18.0        ��b                                                                                                          
       /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/SRC_TUV/ASP//src/EulerianEnvironment/GridPointFields.f90 GRIDPOINTFIELDS       5       ALLOCATEGASCHEMISTRYGRID ALLOCATEGASREACTIONGRID DIFFUSIONCOEFFICIENTOFWATER GASPHASECHEMICALRATES GETAIRDENSITY GETCPM GETDYNAMICVISCOSITYOFAIR GETGASCHEM GETGASBURDEN GETGASCHEMS GETGASCHEMVECTOR GETGASCHEMFROMBURDEN GETM GETMIXINGRATIO GETPARTIALPRESSURE GETPRESS GETPRESSFROMGRID GETRELATIVEHUMIDITY GETRELATIVEHUMIDITYFROMGRID GETRELATIVEHUMIDITYFROMBURDEN GETSATVAPPRESS GETSATMIXINGRATIO GETSATVAPCONCENTRATION GETSATVAPBURDEN GETTEMP GETTEMPFROMGRID GETTHERMALVELOCITYOFAIR GRIDPOINTVOLUME INITIALIZEGRIDPOINTFIELDS MEANFREEPATHOFAIR MOLECULARTHERMALCONDUCTIVITY MOLECULARTHERMALDIFFUSIVITY RESETGASPHASEREACTIONRATES SETGRIDPOINTVOLUMES SETHOMOGENEOUSCHEMFIELDPPB SETHOMOGENEOUSCHEMFIELDMGPERM3 SETHOWMANYGASCHEMS SETWATERVAPORFIELD SETPRESSFIELD SETTEMPFIELD SETYLEN SETHOWMANYEVOLVEGASCHEMS SURFACETENSIONOFWATER SETRELATIVEHUMIDITY SETALLRELATIVEHUMIDITIES UPDATECHEMICALBURDENS UPDATECHEMICALCONCENTRATIONS GETGRIDCELLCHEMBURDEN ADDTOGRIDCELLCHEMBURDEN REPLACEGRIDCELLCHEMBURDEN ALLSAMEGRID GRIDGASCHEM LINEAR2DINTERP #         @                                                       #NUMBCHEMS                                                           #         @                                                       #NUMBRATES                                                           %         @                                                    
                                                                                                 @                                                   
                &                                           %         @                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              %         @                                                    
                                       %         @                              	                     
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     %         @                              
                    
       #SINGLECHEMINDEX              
                                             %         @                                                   
       #SINGLECHEMINDEX              
                                             (        `                                                                      
    p          5 r        5 r                      (        `                                                                      
    p          5 r        5 r                      %         @                                                   
       #BURDEN              
                                      
      %         @                                                   
                               %         @                                                  
       #RHIN                                                               @                                   
       %         @                                                   
       #SINGLECHEMINDEX                                         
                                             %         @                                                   
       %         @                                                    
       %         @                                                   
                                                                       %         @                                                    
                                                                                       %         @                                                   
       #WATERBURDEN                                                                                                                                       
       %         @                                                   
                                %         @                                                     
                                                             %         @                              !                     
                                                                             %         @                               "                     
       %         @                              #                     
       %         @                               $                     
       %         @                              %                     
                                                                                                                                   @                                &     
       #         @                                   '                    #T (   #PRESS )             D @                              (     
                 D @                              )     
       %         @                               *                     
                                                                                                                                  %         @                              +                     
                                                                                     %         @                               ,                     
                                                                                                                       #         @                                   -                    #REACTIONRATES .   #ARRAYSIZE /   #TOTALPEROXY 0   #TOTALACYLPEROXY 1   #H2OCONC 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                 .                   
              &                   &                                                     
                                  /                       p          p            p                                                                    0     
                                                 1     
                                                 2     
       #         @                                  3                     #         @                                   4                    #INDEX 5   #VALUE 6                                                                                                                                 5                                                      6     
       #         @                                   7                    #INDEX 8   #VALUE 9   #MOLECMASS :                                                                                                                                                                                   8                                                      9     
                                                 :     
       #         @                                   ;                    #INHOWMANYGASCHEMS <                                              <            #         @                                   =                    #RH >                                                                                                                >     
       #         @                                  ?                    #PRESS @                                                                              @     
       #         @                                  A                    #TEMP B                                             B     
       #         @                                   C                    #INYLEN D                                              D            #         @                                   E                    #INHOWMANYEVOLVEGASCHEMS F                                              F            %         @                               G                     
                                        #         @                                  H                    #RH I                                                                                                                                               I     
       #         @                                   J                    #RH K             D @                              K     
       #         @                                   L                    #Y M             @                              M                    
     p          5 r N       5 r N                     #         @                                   O                    #Y P                                            P                    
 	    p          5 r N       5 r N                     %         @                               Q                    
       #SINGLECHEMINDEX R             
                                  R           #         @                                   S                    #SINGLECHEMINDEX T   #DELTACHEMBURDEN U             
                                  T                     
                                 U     
      #         @                                  V                    #SINGLECHEMINDEX W   #CHEMBURDEN X             
                                  W                     
                                 X     
                 @                                 Y                     @                                Z                   
                &                                           #         @                                   [                    #X1 \   #Y1 ]   #X2 ^   #Y2 _   #X3 `   #Y3 a             
                                 \     
                
                                 ]     
                
                                 ^     
                
                                 _     
                
                                 `     
                D                                a     
                  @  @                              N                       @  @                                                     @  @                                             �         fn#fn %       b   uapp(GRIDPOINTFIELDS )   6  W       ALLOCATEGASCHEMISTRYGRID 3   �  @   a   ALLOCATEGASCHEMISTRYGRID%NUMBCHEMS (   �  W       ALLOCATEGASREACTIONGRID 2   $  @   a   ALLOCATEGASREACTIONGRID%NUMBRATES ,   d  �       DIFFUSIONCOEFFICIENTOFWATER &     �       GASPHASECHEMICALRATES    �  w      GETAIRDENSITY    
  p       GETCPM )   x
  ^      GETDYNAMICVISCOSITYOFAIR    �  e       GETGASCHEM +   ;  @   a   GETGASCHEM%SINGLECHEMINDEX    {  e       GETGASBURDEN -   �  @   a   GETGASBURDEN%SINGLECHEMINDEX       �       GETGASCHEMS !   �  �       GETGASCHEMVECTOR %   h  \       GETGASCHEMFROMBURDEN ,   �  @   a   GETGASCHEMFROMBURDEN%BURDEN      h       GETM    l  �       GETMIXINGRATIO $   �  @   a   GETMIXINGRATIO%RHIN #   6  �       GETPARTIALPRESSURE 3   �  @   a   GETPARTIALPRESSURE%SINGLECHEMINDEX    �  P       GETPRESS !   F  P       GETPRESSFROMGRID $   �  �       GETRELATIVEHUMIDITY ,   &  �       GETRELATIVEHUMIDITYFROMGRID .   �  �       GETRELATIVEHUMIDITYFROMBURDEN :   {  @   a   GETRELATIVEHUMIDITYFROMBURDEN%WATERBURDEN    �  i       GETSATVAPPRESS "   $  �       GETSATMIXINGRATIO '   �  �       GETSATVAPCONCENTRATION     @  P       GETSATVAPBURDEN    �  P       GETTEMP     �  P       GETTEMPFROMGRID (   0  �       GETTHERMALVELOCITYOFAIR     �  @       GRIDPOINTVOLUME *   1  Z       INITIALIZEGRIDPOINTFIELDS ,   �  @   a   INITIALIZEGRIDPOINTFIELDS%T 0   �  @   a   INITIALIZEGRIDPOINTFIELDS%PRESS "     �       MEANFREEPATHOFAIR -   �  �       MOLECULARTHERMALCONDUCTIVITY ,   t  �       MOLECULARTHERMALDIFFUSIVITY +   4  �      RESETGASPHASEREACTIONRATES 9   	   �   a   RESETGASPHASEREACTIONRATES%REACTIONRATES 5   �   �   a   RESETGASPHASEREACTIONRATES%ARRAYSIZE 7   A!  @   a   RESETGASPHASEREACTIONRATES%TOTALPEROXY ;   �!  @   a   RESETGASPHASEREACTIONRATES%TOTALACYLPEROXY 3   �!  @   a   RESETGASPHASEREACTIONRATES%H2OCONC $   "  H       SETGRIDPOINTVOLUMES +   I"  �       SETHOMOGENEOUSCHEMFIELDPPB 1   �"  @   a   SETHOMOGENEOUSCHEMFIELDPPB%INDEX 1   :#  @   a   SETHOMOGENEOUSCHEMFIELDPPB%VALUE /   z#  �       SETHOMOGENEOUSCHEMFIELDMGPERM3 5   l$  @   a   SETHOMOGENEOUSCHEMFIELDMGPERM3%INDEX 5   �$  @   a   SETHOMOGENEOUSCHEMFIELDMGPERM3%VALUE 9   �$  @   a   SETHOMOGENEOUSCHEMFIELDMGPERM3%MOLECMASS #   ,%  _       SETHOWMANYGASCHEMS 5   �%  @   a   SETHOWMANYGASCHEMS%INHOWMANYGASCHEMS #   �%  �       SETWATERVAPORFIELD &   ^&  @   a   SETWATERVAPORFIELD%RH    �&  t       SETPRESSFIELD $   '  @   a   SETPRESSFIELD%PRESS    R'  R       SETTEMPFIELD "   �'  @   a   SETTEMPFIELD%TEMP    �'  T       SETYLEN    8(  @   a   SETYLEN%INYLEN )   x(  e       SETHOWMANYEVOLVEGASCHEMS A   �(  @   a   SETHOWMANYEVOLVEGASCHEMS%INHOWMANYEVOLVEGASCHEMS &   )  q       SURFACETENSIONOFWATER $   �)  �       SETRELATIVEHUMIDITY '   @*  @   a   SETRELATIVEHUMIDITY%RH )   �*  P       SETALLRELATIVEHUMIDITIES ,   �*  @   a   SETALLRELATIVEHUMIDITIES%RH &   +  O       UPDATECHEMICALBURDENS (   _+  �   a   UPDATECHEMICALBURDENS%Y -   �+  O       UPDATECHEMICALCONCENTRATIONS /   B,  �   a   UPDATECHEMICALCONCENTRATIONS%Y &   �,  e       GETGRIDCELLCHEMBURDEN 6   ;-  @   a   GETGRIDCELLCHEMBURDEN%SINGLECHEMINDEX (   {-  r       ADDTOGRIDCELLCHEMBURDEN 8   �-  @   a   ADDTOGRIDCELLCHEMBURDEN%SINGLECHEMINDEX 8   -.  @   a   ADDTOGRIDCELLCHEMBURDEN%DELTACHEMBURDEN *   m.  m       REPLACEGRIDCELLCHEMBURDEN :   �.  @   a   REPLACEGRIDCELLCHEMBURDEN%SINGLECHEMINDEX 5   /  @   a   REPLACEGRIDCELLCHEMBURDEN%CHEMBURDEN    Z/  @       ALLSAMEGRID    �/  �       GRIDGASCHEM    &0  x       LINEAR2DINTERP "   �0  @   a   LINEAR2DINTERP%X1 "   �0  @   a   LINEAR2DINTERP%Y1 "   1  @   a   LINEAR2DINTERP%X2 "   ^1  @   a   LINEAR2DINTERP%Y2 "   �1  @   a   LINEAR2DINTERP%X3 "   �1  @   a   LINEAR2DINTERP%Y3 &   2  @      HOWMANYEVOLVEGASCHEMS    ^2  @      YLEN     �2  @      HOWMANYGASCHEMS 