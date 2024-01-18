

unset title
set terminal pdf size 10cm,8cm
set output "./out/qnetwork_jonction.pdf"
set grid
set xlabel "time"
set ylabel "discharges"
plot for [j=2:6] "./out/qnetwork_jonction.txt" u ($1-1):(column(j)) w lp lw 2 title "discharges at node".sprintf("%i",j-1)
unset output
unset terminal

unset title
set terminal pdf size 10cm,8cm
set output "./out/vnetwork_jonction.pdf"
set grid
set xlabel "time"
set ylabel "velocities"
plot for [j=2:6] "./out/vnetwork_jonction.txt" u ($1-1):(column(j)) w lp lw 2 title "velocities at node".sprintf("%i",j-1)
unset output
unset terminal




unset title
set terminal pdf size 10cm,8cm
set output "./out/qnetwork_initial.pdf"
set grid
set xlabel "time"
set ylabel "discharges"
plot for [j=2:6] "./out/qnetwork_initial.txt" u ($1-1):(column(j)) w lp lw 2 title "discharges at node".sprintf("%i",j-1)
unset output
unset terminal

unset title
set terminal pdf size 10cm,8cm
set output "./out/vnetwork_initial.pdf"
set grid
set xlabel "time"
set ylabel "velocities"
plot for [j=2:6] "./out/vnetwork_initial.txt" u ($1-1):(column(j)) w lp lw 2 title "velocities at node".sprintf("%i",j-1)
unset output
unset terminal




unset title
set terminal pdf size 10cm,8cm
set output "./out/qnetwork_evaluation.pdf"
set grid
set xlabel "time"
set ylabel "discharges"
plot for [j=2:6] "./out/qnetwork_evaluation.txt" u ($1-1):(column(j)) w lp lw 2 title "discharges at node".sprintf("%i",j-1)
unset output
unset terminal

unset title
set terminal pdf size 10cm,8cm
set output "./out/vnetwork_evaluation.pdf"
set grid
set xlabel "time"
set ylabel "velocities"
plot for [j=2:6] "./out/vnetwork_evaluation.txt" u ($1-1):(column(j)) w lp lw 2 title "velocities at node".sprintf("%i",j-1)
unset output
unset terminal




set terminal pdf size 20cm,15cm font "times,18" enhanced
set output "./out/Tabulated_routing_coefficient_3D_pdf.pdf"
set multiplot layout 2,2
ncol=system("awk 'NR==1{print NF}' ./out/tabulated_gamma_coefficient_3D_pdf.txt")
nrow1=system("awk 'END{print NR}' ./out/tabulated_gamma_coefficient_3D_pdf.txt")
#print nrow
spreading_step=0.1
nblock=26 
set xlabel "Time-steps"
set ylabel "Gamma routing coefficients"
set xrange[0:20]
set xtics 2
set key font "times,14"
do for [j=1:nblock:8]{
    set title "Diffusion=".sprintf("%2.1f",(j*1.-1.)*spreading_step)." s.m^{-1}"
    plot for [i=2:ncol+1:20] "./out/tabulated_gamma_coefficient_3D_pdf.txt" index (j-1) using ($1-1):i w l lw 4 title "Mode = ".sprintf('%2.1f',i*0.1-0.1)
}
unset multiplot
unset output
unset terminal


set terminal pdf size 20cm,15cm font "times,18" enhanced
set output "./out/Routing_coefficient_pdf.pdf"
set multiplot layout 2,2
ncol=system("awk 'NR==1{print NF}' ./out/tabulated_gamma_coefficient_3D_pdf.txt")
set xlabel "Time-steps"
set ylabel "Gamma routing coefficients"
set xrange[0:20]
set xtics 2
set key font "times,14"
spreading_step=0.1
j=25

set title "A) Diffusion=".sprintf("%2.1f",(j*1.-1.)*spreading_step)." s.m^{-1}, dt=900s, dx=1000m"
plot for [i=2:ncol+1:20] "./out/tabulated_gamma_coefficient_3D_pdf.txt" index (j-1) using ($1-1):i w l lw 4 title "Mode = ".sprintf('%2.1f',i*0.1-0.1)
set title "B) Diffusion=".sprintf("%2.1f",(j*1.-1.)*spreading_step)." s.m^{-1}, dt=3600s, dx=1000m"
plot for [i=2:ncol+1:20] "./out/tabulated_gamma_coefficient_3600s_pdf.txt" index (j-1) using ($1-1):i w l lw 4 title "Mode = ".sprintf('%2.1f',i*0.1-0.1)
set title "C) Diffusion=".sprintf("%2.1f",(j*1.-1.)*spreading_step)." s.m^{-1}, dt=900s, dx=250m"
plot for [i=2:ncol+1:20] "./out/tabulated_gamma_coefficient_250m_pdf.txt" index (j-1) using ($1-1):i w l lw 4 title "Mode = ".sprintf('%2.1f',i*0.1-0.1)
set title "D) Diffusion=".sprintf("%2.1f",(j*1.-1.)*spreading_step)." s.m^{-1}, dt=3600s, dx=250m"
plot for [i=2:ncol+1:20] "./out/tabulated_gamma_coefficient_3600s250m_pdf.txt" index (j-1) using ($1-1):i w l lw 4 title "Mode = ".sprintf('%2.1f',i*0.1-0.1)


unset multiplot
unset output
unset terminal






# ~ set terminal pdf size 10cm,8cm
# ~ set output "./out/Gamma_coefficient.pdf"
# ~ set grid
# ~ set xlabel "quantile"
# ~ set ylabel "Gamma coefficients"
# ~ plot "./out/gamma_coefficient_vmax.txt" u ($1-1):2 w lp lw 2 title "vmax", "./out/gamma_coefficient_vmoy.txt" u ($1-1):2 w lp lw 2 title "vmoy", "./out/gamma_coefficient_vmin.txt" u ($1-1):2 w lp lw 2 title "vmin"
# ~ unset output
# ~ unset terminal


# ~ set terminal pdf size 10cm,8cm
# ~ set output "./out/Gamma_coefficient_pdf_cdf.pdf"
# ~ set grid
# ~ set xlabel "quantile"
# ~ set ylabel "Gamma coefficients"
# ~ plot "./out/gamma_coefficient_vmoy.txt" u ($1-1):2 w lp lw 2 title "PDF", "./out/gamma_coefficient_cdf.txt" u ($1-1):2 w lp lw 2 title "CDF"
# ~ unset output
# ~ unset terminal


# ~ set terminal pdf size 10cm,8cm
# ~ set output "./out/Tabulated_routing_coefficient_pdf.pdf"
# ~ ncol=system("awk 'END{print NF}' ./out/tabulated_routing_coefficient_pdf.txt")
# ~ set key font "times,6"
# ~ plot for [i=2:ncol+1:4] "./out/tabulated_routing_coefficient_pdf.txt" u 1:i w l lw 2 title "Mode=".sprintf('%2.1f',i*0.1-0.1)
# ~ unset output
# ~ unset terminal

# ~ set terminal pdf size 10cm,8cm
# ~ set output "./out/New_Tabulated_routing_coefficient_pdf.pdf"
# ~ ncol=system("awk 'END{print NF}' ./out/new_tabulated_routing_coefficient_pdf.txt")
# ~ set key font "times,6"
# ~ plot for [i=2:ncol+1:4] "./out/new_tabulated_routing_coefficient_pdf.txt" u 1:i w l lw 2 title "Mode=".sprintf('%2.1f',i*0.1-0.1)
# ~ unset output
# ~ unset terminal





# ~ set terminal pdf size 10cm,8cm
# ~ set output "./out/Tabulated_routing_coefficient_cdf.pdf"
# ~ ncol=system("awk 'END{print NF}' ./out/tabulated_routing_coefficient_cdf.txt")
# ~ set key font "times,6"
# ~ plot for [i=2:ncol+1:4] "./out/tabulated_routing_coefficient_cdf.txt" u 1:i w l lw 2 title "Mode=".sprintf('%2.1f',i*0.1-0.1)
# ~ unset output
# ~ unset terminal


# ~ set terminal pdf size 10cm,8cm
# ~ set output "./out/Tabulated_Gamma_coefficient_cdf_Igor.pdf"
# ~ ncol=system("awk 'END{print NF}' ./out/tabulated_gamma_coefficient_Igor.txt")
# ~ set key font "times,6"
# ~ plot for [i=2:ncol+1:4] "./out/tabulated_gamma_coefficient_Igor.txt" u 1:i w l lw 2 title "Mode=".sprintf('%2.1f',i*0.1-0.1)
# ~ unset output
# ~ unset terminal

# ~ set terminal pdf size 10cm,8cm
# ~ set output "./out/Tabulated_Gamma_coefficient_cdf.pdf"
# ~ ncol=system("awk 'END{print NF}' ./out/tabulated_gamma_coefficient_cdf.txt")
# ~ set key font "times,6"
# ~ plot for [i=2:ncol+1:4] "./out/tabulated_gamma_coefficient_cdf.txt" u 1:i w l lw 2 title "Mode=".sprintf('%2.1f',i*0.1-0.1)
# ~ unset output
# ~ unset terminal

# ~ set terminal pdf size 10cm,8cm
# ~ set output "./out/Tabulated_Gamma_coefficient_pdf.pdf"
# ~ ncol=system("awk 'END{print NF}' ./out/tabulated_gamma_coefficient_pdf.txt")
# ~ set key font "times,6"
# ~ plot for [i=2:ncol+1:4] "./out/tabulated_gamma_coefficient_pdf.txt" u 1:i w l lw 2 title "Mode=".sprintf('%2.1f',i*0.1-0.1)
# ~ unset output
# ~ unset terminal

# ~ set terminal pdf size 10cm,8cm
# ~ set output "./out/Interpolated_Gamma_coefficient.pdf"
# ~ set grid
# ~ set xlabel "quantile"
# ~ set ylabel "Gamma coefficients"
# ~ plot "./out/gamma_coefficient_true.txt" u ($1-1):2 w lp lw 2 title "True", "./out/interpolated_gamma_coefficient.txt" u ($1-1):2 w lp lw 2 title "Interpolated"
# ~ unset output
# ~ unset terminal

# ~ unset title
# ~ set terminal pdf size 10cm,8cm
# ~ set output "./out/qnetwork.pdf"
# ~ set grid
# ~ set xlabel "time"
# ~ set ylabel "discharges"
# ~ plot for [j=2:6] "./out/qnetwork.txt" u ($1-1):(column(j)) w lp lw 2 title "discharges at node".sprintf("%i",j-1)
# ~ unset output
# ~ unset terminal

# ~ unset title
# ~ set terminal pdf size 10cm,8cm
# ~ set output "./out/vnetwork.pdf"
# ~ set grid
# ~ set xlabel "time"
# ~ set ylabel "velocities"
# ~ plot for [j=2:6] "./out/vnetwork.txt" u ($1-1):(column(j)) w lp lw 2 title "velocities at node".sprintf("%i",j-1)
# ~ unset output
# ~ unset terminal


