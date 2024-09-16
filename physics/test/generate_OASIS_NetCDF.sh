#!/usr/bin/bash

InVars=(I_SST:-1.84 I_SSS:28 I_Uocn:-0.14 I_Vocn:0.71)
# Optionally add
# InVars+=(I_SSH:14.8 I_MLD:14.8)

OutVars=(I_taux I_tauy I_taumod I_fwflux I_rsnos I_rsso I_sfi I_conc)

cat > namcouple <<EOF
#
# OASIS Input Output basic test for NextSIM-DG
#
#########################################################################
\$NFIELDS
 $(( ${#InVars[@]} + ${#OutVars[@]} ))
\$END
###########################################################################
\$RUNTIME
# The total simulated time for this run in seconds
  1
\$END
###########################################################################
\$NLOGPRT
# Amount of information written to OASIS3-MCT log files (see User Guide)
  0
\$END
###########################################################################
\$STRINGS
# updateBefore fields read in from file
EOF

for var in ${InVars[@]} ; do
    av=(${var//:/ })
    field=${av[0]}
    value=${av[1]}

    cat > tmp.cdl <<EOF
netcdf ${field}_In {
dimensions:
        nx = 1 ;
        ny = 1 ;
        time = UNLIMITED ; // (1 currently)
variables:
        int time(time) ;
        double ${field}(time, ny, nx) ;
data:

 time = 0 ;

 ${field} = ${value} ;
}
EOF

    ncgen -b -o ${field}_In.nc tmp.cdl
    rm -f tmp.cdl

echo "   ${field} ${field} 1 1 0 ${field}_In.nc INPUT" >> namcouple

done

cat >> namcouple <<EOF

# updateAfter field written to file
EOF

for var in ${OutVars[@]} ; do
cat >> namcouple <<EOF
   ${var} ${var} 1 1 0 null OUTPUT
     1 1 1 1 null null
EOF


done

echo "\$END" >> namcouple
