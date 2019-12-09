#Importing packages
import os, time
from datetime import datetime
from shutil import rmtree

#######################
#MODIFIED  BY the USER#
################################################################
Short_name_station = "DDU" #Define folter station name
name_station = '_'.join(Short_name_station.split(' '))
path = "Data/"+name_station+"/"
year =  2017
month = 2
day  = 9

#Send email

NC2transfer1 = "Daily-schedules/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)+"_sub.nc"
NC2transfer2 = "Daily-schedules/"+name_station+"_"+str(year)+str(month).zfill(2)+str(day).zfill(2)+".png"
                
sendEmail = '/'.join(str(os.path.abspath('lib/Sendemail/sendEmail.exe')).split('\\'))

attached1 = '/'.join(str(os.path.abspath(NC2transfer1)).split('\\'))
attached2 = '/'.join(str(os.path.abspath(NC2transfer2)).split('\\'))

FROM = "apres3@ddu.ipev.fr"
TO = "claudio.duran@univ-grenoble-alpes.fr"
SERVER = "smtp.ddu.ipev.fr:25"
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y")
MMS = "Daily MRR report from DDU, Antarctica, for the past UTC day "+dt_string

os.system(sendEmail+' -a '+attached1+' '+attached2+' -f '+FROM+' -t '+TO+' -u "DAILY MRR REPORT" -s '+SERVER+' -m '+MMS)
time.sleep(15)
os.remove(attached1)
os.remove(attached2)
rmtree("Data/DDU/MK_processed")
rmtree("Data/DDU/temp")
rmtree("Data/DDU/Plots")

