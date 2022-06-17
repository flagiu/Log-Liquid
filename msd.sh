make msd
./msd-log.x
if(($?!=0)); then exit 1; fi
./msd-lin.x
if(($?!=0)); then exit 1; fi
for NAT in {1,3}
do
	cat msd_nat${NAT}-log.dat msd_nat${NAT}-lin.dat >> /dev/null 2>> /dev/null
	if(($?==0)) #check if those file were created
	then
		cat msd_nat${NAT}-log.dat msd_nat${NAT}-lin.dat > msd_nat${NAT}.dat
		cat msdcm_nat${NAT}-log.dat msdcm_nat${NAT}-lin.dat > msdcm_nat${NAT}.dat
		paste msd_nat${NAT}.dat msdcm_nat${NAT}.dat | awk '{printf "%f %f %f\n",$1,$2-$5,sqrt($3**2+$6**2)}' > msdsub_nat${NAT}.dat
		echo "Created msdsub_nat${NAT}.dat"
	fi
done
echo ""
