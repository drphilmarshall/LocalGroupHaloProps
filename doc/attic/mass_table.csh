# Make latex table of masses:

#  Columns: M_LG, M_MW, M_M31, M_M33, A200, M_M31/M_MW 
#  Rows:    
#  1  isolated pairs
#  2  + M31 D and vr
#  3  + MW mass
#    -----------------
#  4  isolated triplets
#  5  + M31 D and vr
#  6  + M33 D and vr
#  7  + MW mass

set klobber = 0

set files = (\
LG_pairs_prior.cpt \
LG_pairs_M31_D+vr.cpt \
LG_pairs_M31_D+vr_MW-Mvir.cpt \
LG_triplets_prior.cpt \
LG_triplets_M31_D+vr.cpt \
LG_triplets_M33+M31_D+vr.cpt \
LG_triplets_M33+M31_D+vr_MW-Mvir.cpt )

# fig2c_pairs_masses.log \
# fig2c_triplets_masses.log \
# fig2d_pairs_masses.log \
# fig2d_triplets_masses.log \
# fig3a_pairs_mass_comparison.log \
# fig3a_triplets_mass_comparison.log )

set triplets = ( \
0 \
0 \
0 \
1 \
1 \
1 \
1 )

set info = ( \
'Pairs              ' \
'$+$ M31 $D,v$      ' \
'$+$ Bolshoi MW mass' \
'Triplets           ' \
'$+$ M31 $D,v$      ' \
'$+$ M33 $D,v$      ' \
'$+$ Bolshoi MW mass' \
)

# Table header:

set tabfile = ~/Dropbox/LocalGroup/doc/mass_table.tex
echo "Writing table to $tabfile"

echo \
'\begin{deluxetable*}{lccccccc}\
\tabletypesize{\small}\
\tablecaption{\label{tab:masses}\
Local group mass estimates}\
\tabletypesize{\small}\
\tablehead{\
Constraints & \
$\log_{10} M_{\rm MW} / M_{\odot}$       & \
$\log_{10} M_{\rm M31} / M_{\odot}$      & \
$\log_{10} M_{\rm M31} / M_{\rm MW}$     & \
$\log_{10} M_{\rm M33} / M_{\rm MW}$     & \
$\log_{10} M\prime_{\rm LG} / M_{\odot}$ & \
$\log_{10} M_{\rm TA} / M_{\odot}$       & \
$\log_{10} A_{200}$ \\} \
\startdata' > $tabfile


# Pairs:

\rm -f junk
foreach i ( `seq $#files` )

  set cptfile = $files[$i]

# Compute estimates and store in a file:

  set estfile = $cptfile:r.est
  if ($klobber) \rm -f $estfile
  if (! -e $estfile) then
    PointEstimator.py -w 1 $cptfile > $estfile
    wc -l $estfile
  endif  

# Extract estimates and format them in the table:
#   Columns: M_LG, M_MW, M_M31, M_M33, M_TA, A200, M_M31/M_MW 

  set MLG = `grep -e M_ $estfile | grep -e LG | grep -e prime | cut -d'=' -f2 | sed s/'{'/%/g | sed s/'}'/=/g`
  set MMW = `grep -e M_ $estfile | grep -e MW | grep -e odot | cut -d'=' -f2 | sed s/'{'/%/g | sed s/'}'/=/g` 
  set M31 = `grep -e M_ $estfile | grep -e M31 | grep -e odot | cut -d'=' -f2 | sed s/'{'/%/g | sed s/'}'/=/g` 
  set ratio = `grep -e M_ $estfile | grep -e M31 | grep -e MW | cut -d'=' -f2 | sed s/'{'/%/g | sed s/'}'/=/g`
  
  if $triplets[$i] then
    set M33 = `grep -e M_ $estfile | grep -e M33 | grep -e odot | cut -d'=' -f2 | sed s/'{'/%/g | sed s/'}'/=/g` 
  else 
    set M33 = '--'
  endif
  
  set MTA = `grep -e M_ $estfile | grep -e TA | grep -e odot | cut -d'=' -f2 | sed s/'{'/%/g | sed s/'}'/=/g`
  set A200 = `grep -e A_ $estfile | grep -e 200 | cut -d'=' -f2 | sed s/'{'/%/g | sed s/'}'/=/g`

  echo "$info[$i]  &  $MMW  &  $M31  &  $ratio  &  $M33  &  $MLG  &  $MTA  &   $A200 \\" >> junk
 
end  

cat junk | sed s/=/'}'/g | sed s/%/'{'/g >> $tabfile
\rm -f junk

echo '\enddata\
\tablecomments{}\
\end{deluxetable*}' >> $tabfile
