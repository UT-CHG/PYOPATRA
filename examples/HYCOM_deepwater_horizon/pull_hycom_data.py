from ftplib import FTP
from datetime import date, timedelta
from time import time

if __name__ == '__main__':
    total_days = 8 * 7
    start_date = date(2010, 4, 20)
    end_date = start_date + timedelta(days=total_days)

    print("Fetching HYCOM data from {} to {}".format(start_date, end_date))

    times = ['000', '003', '006', '009', '012', '015', '018', '021']

    with FTP('ftp.hycom.org') as ftp:
        ftp.login()
        ftp.cwd('datasets/GOMu0.04/expt_50.1/data/netcdf/2010/')

        for day_since_start in range(total_days):
            date = start_date + timedelta(days=day_since_start)

            for time_str in times:
                file = 'hycom_gomu_501_{}{:02d}{:02d}00_t{}.nc'.format(date.year, date.month, date.day, time_str)

                start = time()
                with open('data/' + file, 'wb') as fp:
                    ftp.retrbinary('RETR {}'.format(file), fp.write)

                print('Downloaded {} in {} seconds'.format(file, time() - start))
