#!/usr/bin/env bash

file="$1"
#f2008
keywords=(block codimension concurrent contiguous critical error stop submodule sync all sync images sync memory lock unlock)
#f2003
keywords+=(abstract associate asynchronous bind class deferred enum enumerator extends final flush generic import non_overridable nopass pass protected value volatile wait)
#f95
keywords+=(elemental forall pure)
#f90
keywords+=(allocatable allocate case contains cycle deallocate elsewhere exit include interface intent module namelist nullify only operator optional pointer private procedure public recursive result select sequence target use while where)
#f77
keywords+=(assign backspace block data call close common continue data dimension do else end endfile endif entry equivalence external format function goto if implicit inquire intrinsic open parameter pause print program read return rewind rewrite save stop subroutine then write)

keywords+=(none impure len kind)
#noy keywords but I am too lazy
keywords+=(intrinsic external none procedure c_ptr c_funptr asynchronous nopass non_overridable pass protected volatile extends
import non_intrinsic value bind deferred generic final enumerator contiguous)
keywords+=(access blank direct exist file fmt form formatted iostat name named nextrec number opened rec recl sequential status
unformatted unit format pad position action delim readwrite eor advance nml namelist decimal round iomsg newunit)
keywords+=(alog alog10 amax0 amax1 amin0 amin1 amod cabs ccos cexp clog csin csqrt dabs dacos dasin datan datan2 dcos dcosh ddim
dexp dint dlog dlog10 dmax1 dmin1 dmod dnint dsign dsin dsinh dsqrt dtan dtanh float iabs idim idint idnint ifix isign max0 max1
min0 min1 sngl abs acos aimag aint anint asin atan atan2 char cmplx conjg cos cosh exp ichar index int log log10 max min nint sign
sin sinh sqrt tan tanh dim lge lgt lle llt mod adjustl adjustr all allocated any associated bit_size btest ceiling count cshift
date_and_time digits dot_product eoshift epsilon exponent floor fraction huge iand ibclr ibits ibset ieor ior ishft ishftc lbound
len_trim matmul maxexponent maxloc maxval merge minexponent minloc minval modulo mvbits nearest pack precision present product radix
random_number random_seed range repeat reshape rrspacing scale scan selected_int_kind selected_real_kind set_exponent shape size
spacing spread sum system_clock tiny transpose trim ubound unpack verify achar iachar transfer dble dprod null cpu_time
command_argument_count get_command get_command_argument get_environment_variable is_iostat_end is_iostat_eor move_alloc new_line
selected_char_kind same_type_as extends_type_of iso_c_binding c_loc c_funloc c_associated c_f_pointer c_f_procpointer
ieee_arithmetic ieee_support_underflow_control ieee_get_underflow_mode ieee_set_underflow_mode acosh asinh atanh bessel_j0 bessel_j1
bessel_jn bessel_y0 bessel_y1 bessel_yn erf erfc erfc_scaled gamma log_gamma hypot norm2 atomic_define atomic_ref
execute_command_line leadz trailz storage_size merge_bits bge bgt ble blt dshiftl dshiftr findloc iall iany iparity image_index
lcobound ucobound maskl maskr num_images parity popcnt poppar shifta shiftl shiftr this_image)

keywords+=(character complex integer character real logical type in out elemental pure impure abstract class associate enum)

special=(interface module procedure program associate select block submodule subroutine critical type do where enum file function forall if)

special2=("else where" "else if" "in out" "go to" "select case" "select type")

for i in "${special2[@]}"; do
  sed -i -e "s/\b\(${i/ /}\)\b/$i/I" $file
done 

for i in ${special[@]}; do
  sed -i -e "s/\b\(end${i}\)\b/end $i/I" $file
done 
for i in ${keywords[@]}; do 
  sed -i -e "s/\b\(${i}\)\b/\L\u\1/I" $file
done
sed -i -e "s/\b\(In Out\)\b/InOut/I" $file
