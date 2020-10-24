#!/usr/bin/env python3
import argparse
import xlrd
from pathlib import Path
import matplotlib as mpl

mpl.use('Agg')  # standard rendering tool for matplotlib above
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from scipy import optimize

ERROR_LOG = Path("errors.log")


def log_error(message):
    with ERROR_LOG.open("a") as log:
        log.write(str(message))


def main():
    # parse args
    parser = argparse.ArgumentParser(description="Run melatonin analysis.")
    parser.add_argument("dataset", type=str, help="path to the dataset")
    parser.add_argument("threshold", type=float, help="threshold value (10 <= x <= 90)")
    args = parser.parse_args()

    # validate args
    if args.dataset is None:
        raise Exception("No dataset specified.")
    if args.threshold is None:
        raise Exception("No threshold specified.")
    dataset = Path(args.dataset)
    if not dataset.is_file():
        raise Exception(f"Dataset not a valid file: {dataset}")
    threshold = args.threshold
    if not 10 <= threshold <= 90:
        raise Exception(f"Invalid threshold: {threshold}")
    dataset_name = dataset.stem
    destpath = dataset.parent.joinpath(get_output_filename(dataset_name, threshold))

    # parse the dataset
    time, data = parse_workbook(xlrd.open_workbook(str(dataset)))

    # run calculation and generate a plotted image
    run(dataset_name, destpath, time, data, threshold)


def get_output_filename(dataset_name, threshold):
    return f"{dataset_name}_{threshold:,.1f}.png"


def parse_workbook(workbook):
    """
    Parse an Excel spreadsheet for time, data.

    :param workbook:    The spreadsheet object.
    :return:            (time, data)
    """
    sheet = workbook.sheet_by_index(0)

    # retrieve time (col A) and data (col B)
    time = []
    data = []
    for row in range(sheet.nrows):
        time.append(sheet.cell_value(row, 0))
        data.append(sheet.cell_value(row, 1))

    return time, data


def cos_fit(params, x):
    """
    A cosine curve fitting function.

    :param params: The function parameters.
    :param x:      The time array.
    :return:       The solution.
    """
    h, b, v, p = params
    return h * (np.cos((2 * np.pi * (x + v)) / p)) + b


def residuals(params, x, y):
    """
    The curve fitting residuals.

    :param params:  The curve function parameters.
    :param x:       The time array.
    :param y:       The data array.
    :return:        The residuals.
    """
    return cos_fit(params, x) - y


def nonlinear_regression(x, y, params_guess, bounds, max_nfev):
    """
    Calculate nonlinear regression via least-squares.  Raise an error if curve fitting failed.

    :param x:            The time array.
    :param y:            The data array.
    :param params_guess: The initial params guess.
    :param bounds:       The solution bounds.
    :param max_nfev:     The maximum number of functional evaluations before giving up.
    :return:             The solved params and residuals.
    """
    result = optimize.least_squares(residuals, params_guess, loss="linear", bounds=bounds, max_nfev=max_nfev,
                                    args=(x, y))
    if not result.success:
        raise Exception("Curve fit error: " + result.message)
    return result.x, result.fun


def pearson(x, y):
    """
    Computes pearson correlation coefficient between arrays x and y.

    :param x:   The fitted curve data.
    :param y:   The measured data.
    :return:    The pearson correlation coefficient.
    """
    xlen = len(x)
    ylen = len(y)
    xsum = x.sum(dtype=float)
    ysum = y.sum(dtype=float)

    # means
    xmean = float(xsum / xlen)
    ymean = float(ysum / ylen)

    # sum of squares
    ss_x = float(sum((i - xmean) ** 2 for i in x))
    ss_y = float(sum((i - ymean) ** 2 for i in y))

    # standard deviation
    sd_x = float((ss_x / xlen) ** .5)
    sd_y = float((ss_y / ylen) ** .5)

    # covariance
    cov_xy = float(sum((a - xmean) * (b - ymean) for a, b in zip(x, y)))

    # check for divide-by-zero
    if sd_x == 0:
        raise Exception("Unable to calculate the Pearson correlation coefficient without dividing by zero (sd_x = 0)")
    if sd_y == 0:
        raise Exception("Unable to calculate the Pearson correlation coefficient without dividing by zero (sd_y = 0)")
    return float((cov_xy / xlen) / (sd_x * sd_y))


def peak_time_data(all_time, all_data, index_coords, onset_coords):
    """
    Filters all_time, and all_data arrays into arrays of lists, each sub-list beginning at onset, ending at offset
    all other coordinates ignored.
    Time list written to peak_time_list, data written to peak_data_list.

    :param all_time:
    :param all_data:
    :param index_coords:
    :param onset_coords:
    :return:
    """
    cos_peak_data_list = []
    cos_peak_time_list = []
    step = 0
    while step < (len(index_coords)):
        if len(index_coords) <= 1:
            log_error("Only 1 mesor crossing; cannot compute peak duration")
            step = len(index_coords)
        elif index_coords[0] == onset_coords[0] and len(index_coords) % 2 == 0 and len(index_coords) >= 2:
            cos_peak_data_list.append(all_data[index_coords[step]:index_coords[step + 1] + 1])
            cos_peak_time_list.append(all_time[index_coords[step]:index_coords[step + 1] + 1])
            step += 2
        elif index_coords[0] == onset_coords[0] and len(index_coords) % 2 != 0 and len(index_coords) >= 3:
            index_coords.pop()
            cos_peak_data_list.append(all_data[index_coords[step]:index_coords[step + 1] + 1])
            cos_peak_time_list.append(all_time[index_coords[step]:index_coords[step + 1] + 1])
            step += 2
        elif index_coords[0] != onset_coords[0] and len(index_coords) % 2 == 0 and len(index_coords) >= 3:
            index_coords = index_coords[1:-1]
        elif index_coords[0] != onset_coords[0] and len(index_coords) % 2 != 0 and len(index_coords) >= 3:
            index_coords = index_coords[1:]
        elif index_coords[0] != onset_coords[0] and len(index_coords) % 2 == 0 and len(index_coords) <= 2:
            raise Exception("Not enough coordinates to determine")
        else:
            print("On-off coordinates not found")
            step = len(index_coords)
    return cos_peak_time_list, cos_peak_data_list


def midpoint_peak_auc(time, data):
    """
    Midpoint auc calculation.

    :param time:    The time array.
    :param data:    The data array.
    :return:
    """
    total_sum = 0
    for i in range(len(time) - 1):
        mp_sum = (time[i + 1] - time[i]) * ((data[i] + data[i + 1]) / 2)
        total_sum += mp_sum
    return total_sum


def get_real_cycles(all_time, onset_coords, offset_coords, fitted_on, fitted_off, day_cycles, p):
    """
    Get the onset/offset crossings of the data curve which correlate to those of the fitted curve.

    :param all_time:        The data's time coordinates, plus those of its curve's crossings.
    :param onset_coords:    The indexes of onset crossings.
    :param offset_coords:   The indexes of offset crossings.
    :param fitted_on:       The time coordinate of the first onset crossing of the fitted curve.
    :param fitted_off:      The time coordinate of the first offset crossing of the fitted curve.
    :param day_cycles:      The number of day-cycles.
    :param p:               The period of the fitted curve.
    :return:                (real_onset_coords, real_offset_coords)
    """
    slightly_less_than_half_of_fitted_duration = abs((fitted_off - fitted_on) / 2.2)

    # generate fitted onsets/offsets, extrapolating one additional cycle in both directions to cover edge cases
    fitted_onsets = [fitted_on + p * day for day in range(-1, day_cycles + 1)]
    fitted_offsets = [fitted_off + p * day for day in range(-1, day_cycles + 1)]

    # create an ordered indexing of all fitted crossings
    all_fitted_xings = sorted(fitted_onsets + fitted_offsets)
    all_fitted_xings_indexed = {i: all_fitted_xings[i] for i in range(len(all_fitted_xings))}

    def find_real_coords(xing_coords, fitted_xings, min_or_max):
        # find the closest fitted xing for each xing
        xing_map = {}
        for xing_coord in xing_coords:
            xing = all_time[xing_coord]
            distances_to_fitted_xings = {fitted_xing: abs(xing - fitted_xing) for fitted_xing in fitted_xings}
            xing_map[xing_coord] = min(distances_to_fitted_xings, key=distances_to_fitted_xings.get)

        # for each fitted xing:
        #   get a list of xings which claim it as the closest one,
        #   run those claimants through a validity filter,
        #   find the best claimant using the given min or max function
        best_claimants = []
        claimant_groups = {}
        for fitted_xing in fitted_xings:
            fitted_xing_index = all_fitted_xings.index(fitted_xing)
            prev_fitted_xing = all_fitted_xings_indexed.get(fitted_xing_index - 1)
            next_fitted_xing = all_fitted_xings_indexed.get(fitted_xing_index + 1)

            claimants = {xing_coord: all_time[xing_coord] for xing_coord in xing_map
                         # is this the closest one?
                         if fitted_xing == xing_map[xing_coord]
                         # validity filter
                         and is_valid_claimant(xing_coord, prev_fitted_xing, next_fitted_xing)}
            if claimants:
                best_claimant = min_or_max(claimants, key=claimants.get)
                best_claimants.append(best_claimant)
                claimant_groups[best_claimant] = (sorted(claimants.keys()))
        return best_claimants, claimant_groups

    def is_valid_claimant(xing_coord, prev_fitted_xing, next_fitted_xing):
        # reject claimants which are beyond the previous or next fitted crossings,
        # or within `fitted_duration / 2` of those crossings
        valid = True
        xing = all_time[xing_coord]
        if prev_fitted_xing is not None:
            if xing < (prev_fitted_xing + slightly_less_than_half_of_fitted_duration):
                valid = False
        if next_fitted_xing is not None:
            if xing > (next_fitted_xing - slightly_less_than_half_of_fitted_duration):
                valid = False
        return valid

    def remove_sandwiches(xing_coords, first_slice, last_slice):
        # remove xing dips that are sandwiched within claimant groups
        sandwiches = [xing_coord for xing_coord in xing_coords
                      if first_slice < xing_coord < last_slice]
        for sandwich in sandwiches:
            xing_coords.remove(sandwich)

    # use the earliest onset (min)
    real_onset_coords, onset_claimant_groups = find_real_coords(onset_coords, fitted_onsets, min)

    # use the latest offset (max)
    real_offset_coords, offset_claimant_groups = find_real_coords(offset_coords, fitted_offsets, max)

    # filter out offsets which are between a valid onset and its highest fellow claimant
    for real_onset_coord in real_onset_coords:
        onset_claimant_group = onset_claimant_groups[real_onset_coord]
        first_onset = real_onset_coord
        last_onset = onset_claimant_group[-1]
        remove_sandwiches(real_offset_coords, first_onset, last_onset)

    # filter out onsets which are between a valid offset and its lowest fellow claimant
    for real_offset_coord in real_offset_coords:
        offset_claimant_group = offset_claimant_groups[real_offset_coord]
        first_offset = offset_claimant_group[0]
        last_offset = real_offset_coord
        remove_sandwiches(real_onset_coords, first_offset, last_offset)

    # for single-days, check for omissions due to curve shift
    if day_cycles == 1:
        if len(real_onset_coords) == 0 and len(onset_coords) > 0:
            # use the first onset if it is less-than the first offset
            if len(offset_coords) == 0 or onset_coords[0] < min(offset_coords):
                real_onset_coords.append(onset_coords[0])
        if len(real_offset_coords) == 0 and len(offset_coords) > 0:
            # use the last offset if it is greater-than the last onset
            if len(onset_coords) == 0 or offset_coords[-1] > max(onset_coords):
                real_offset_coords.append(offset_coords[0])

    return real_onset_coords, real_offset_coords


def run(name_of_dataset, destpath, x, y, ask_user_threshold):
    # wrap the measured data in numpy arrays
    time = np.array(x)
    data = np.array(y)

    # init curve fitting params
    initial_params_guess = [np.max(data), np.min(data), -3, 24]
    max_nfev = 10000000000
    bounds = ([-np.inf, np.min(data), -np.inf, 3], [np.inf, np.inf, np.inf, np.inf])

    # find fitted curve
    cos_params, res = nonlinear_regression(time, data, initial_params_guess, bounds, max_nfev)

    # calculate the pearson correlation coefficient (r) and the coefficient of determination (r^2)
    r = pearson(cos_fit(cos_params, time), data)
    r2 = r ** 2

    # calculate the sum of squared residuals (ss)
    # the defining value of the "fitted" function is to return the smallest possible ss
    ss = sum((cos_fit(cos_params, time) - data) ** 2)
    df = len(data) - 2
    RSD = np.sqrt(ss / df)

    # calculate acrophase
    # (x, y coordinates of highest point of cosine function: lsq_peak value is y value, acro(x) finds x
    cos_peak_value = cos_params[1] + cos_params[0]

    # acrophase time directly related to cos_params[2], "v" which is x offset for fit.
    # for example if v = 0, peak timing would occur at x = 0
    def acro(x):
        if x < 0:
            return abs(x) + cos_params[3]
        else:
            return cos_params[3] - x

    # find mesor (midpoint between peak and trough, also when standard cosine function cos(x) for x = -pi/2 and x = pi/2
    # but most easily is found from cos_params[1], "b" as the vertical offset raising or lowering ht of function)

    # determine threshold where 0 = min (trough) value, 50 = mesor, 100 = max (peak) value
    user_threshold = cos_params[1]
    cos_range = 2 * (cos_peak_value - cos_params[1])
    lsq_min = cos_peak_value - cos_range
    if 0 <= ask_user_threshold <= 100:
        user_threshold = ((ask_user_threshold / 100) * cos_range) + lsq_min
    else:
        user_threshold = cos_params[1]
        print("Threshold out of range. Default of 50 (100% mesor) used.")

    # create list populated from mesor values equal in length to time (and thus data) list for later comparison
    cos_mesor = user_threshold
    cos_mesor_list = [cos_mesor] * len(time)

    # acrophase to line up with appropriate time series (ie if dataset time starts after hr 30+)
    cos_acro = acro(cos_params[2])
    while cos_acro <= np.min(time):
        cos_acro += cos_params[3]

    # cos_std_cos removes cos_params[1] "b" and cos_params[0] "h" to return y values back to -1<=y<=1
    cos_std_cos = (user_threshold - cos_params[1]) / cos_params[0]

    # cos_cos_prime returns the first order derivatives to cos_std_cos function
    cos_cos_prime = -np.sin(cos_std_cos)

    # cos_acos returns the cos integral of cos_std_cos
    # assumed beneficial for known y values where it is desired to know x (time)
    cos_acos = np.arccos(cos_std_cos)

    # cos_time_up generates timepoints where mesor is crossed by cos_fit going up (onset)
    cos_time_up = (((2 * np.pi - cos_acos) * cos_params[3]) / (2 * np.pi)) - cos_params[2]

    # cos_time_down generates timepoints where mesor is crossed by cos_fit going down (offset)
    cos_time_down = ((cos_acos * cos_params[3]) / (2 * np.pi)) - cos_params[2] + cos_params[3]

    cos_fitted_on = cos_time_up
    cos_fitted_off = cos_time_down
    cos_fitted_duration = cos_fitted_off - cos_fitted_on
    cos_fitted_on = cos_acro - .5 * cos_fitted_duration
    cos_fitted_off = cos_acro + .5 * cos_fitted_duration
    # tuple pairing time and mesor values
    cos_mesor_tuple = list(zip(time, cos_mesor_list))

    # tuple pairing time and data values
    cos_time_data_tuple = list(zip(time, data))

    # repeat acro coords across figure rather than displaying 1 point \
    # should plot point for every peak occurring within dataset
    cos_acro_list_x = []
    cos_acro_list_y = []
    cos_acro_point_total = cos_acro
    while cos_acro_point_total <= np.min(time):
        cos_acro_point_total += cos_params[3]

    while cos_acro_point_total <= np.max(time):
        cos_acro_list_x.append(cos_acro_point_total)
        cos_acro_list_y.append(cos_peak_value)
        cos_acro_point_total += cos_params[3]

    # mesor-data intersection points
    # below returns data index points prior to mesor crossing, writes these index values to cos_y_int
    idx = np.argwhere(np.diff(np.sign(cos_mesor_list - data))).flatten()
    cos_x_int = time[idx]
    cos_y_int = [cos_mesor] * len(cos_x_int)

    # finds time values where mesor intersects with data (assuming straight line from point to point)
    # by generating straight line y = mx+b from two known points per interval, and inverse to solve for y w/known time
    # writes each intersection timepoint to cos_crossing_points list (y value is always mesor, this list is x "time")
    cos_crossing_points = []
    for intersect in idx:
        cos_crossing_points.append(
            (cos_mesor - ((data[intersect]) - ((data[intersect + 1] - data[intersect]) / (time[intersect + 1] -
                                                                                          time[intersect])) *
                          time[intersect])) / ((data[intersect + 1] - data[intersect]) / (time[intersect + 1] -
                                                                                          time[intersect])))

    cos_crossing_points = np.array(cos_crossing_points)
    cos_crossing_points_unstring = ", ".join("{0:.3f}".format(num) for num in cos_crossing_points)

    # determine area under curve between onset and offset using trapezoidal rule
    # combine time,data coords with cos_crossing_points, cos_mesor list, then arrange by time
    cos_crossing_coords = sorted(list(zip(cos_crossing_points, cos_mesor_list)) + cos_time_data_tuple)
    cos_sorted_coords = list(zip(*cos_crossing_coords))
    # cos_all_time is every original timepoint plus mesor intersection timepoints
    cos_all_time = cos_sorted_coords[0]
    # cos_all_data is every original datapoint plus mesor value when mesor intersects data
    cos_all_data = cos_sorted_coords[1]

    # find position in sorted coords where crossing points appear (start and end of each auc computation)
    cos_index_coords = []
    for pts in cos_crossing_points:
        for items in cos_all_time:
            if items == pts:
                cos_index_coords.append(cos_all_time.index(items))

    # introduce lists that will be written as a function of whether curve is going up during mesor crossing (onset)
    # or down during crossing (offset) These are cos_onset_coords and cos_offset_coords, respectively.
    # onset and offset cos_index_coords display the index number of these locations relative to position in cos_all_time
    cos_onset_coords = []
    cos_onset_cos_index_coords = []
    cos_offset_coords = []
    cos_offset_cos_index_coords = []
    for pts in cos_index_coords:
        if cos_all_data[pts + 1] > cos_all_data[pts]:
            cos_onset_coords.append(pts)
            cos_onset_cos_index_coords.append(cos_index_coords.index(pts))

    # for pts in cos_index_coords[1:len(cos_index_coords)]:
    for pts in cos_index_coords:
        if cos_all_data[pts + 1] < cos_all_data[pts]:
            cos_offset_coords.append(pts)
            cos_offset_cos_index_coords.append(cos_index_coords.index(pts))

    # find beginning and end of first period interval (cycle 1)
    cos_day_1_interval = []
    cos_day_1_interval.append(cos_acro - .5 * cos_params[3])
    cos_day_1_interval.append(cos_acro + .5 * cos_params[3])

    # find number of cycles in dataset
    cos_day_cycles = 1
    cos_day_step = cos_params[3]
    cos_day_max = cos_acro + .5 * cos_params[3]
    while cos_day_max < np.max(time):
        cos_day_cycles += 1
        cos_day_max += cos_day_step
        cos_day_1_interval.append(cos_day_max)

    # breakdown number of time & data values per cycle, with cycle index
    cos_cycle_time = []
    cos_cycle_data = []
    cos_cycle_index = []
    cos_cycle_time_index = []
    cos_cycle_step = 1
    while cos_cycle_step < len(cos_day_1_interval) - 1:
        for hr in time:
            if cos_day_1_interval[cos_cycle_step - 1] <= hr <= cos_day_1_interval[cos_cycle_step]:
                cos_cycle_time.append(hr)
                cos_cycle_index.append(cos_cycle_step)
                cos_cycle_time_index.append(x.index(hr))
            elif hr > cos_day_1_interval[cos_cycle_step]:
                cos_cycle_step += 1
                cos_cycle_time.append(hr)
                cos_cycle_index.append(cos_cycle_step)
                cos_cycle_time_index.append(x.index(hr))

    cos_sample_time_per_day_list = []
    if cos_day_cycles == 1:
        cos_sample_time_per_day_list.append(len(time))
        for hr in time:
            cos_cycle_index.append(1)

    # generate list containing number of datapoints per cycle
    cos_cycle_step = 1
    cos_cycle_time_holder = []
    cos_cycle_data_holder = []
    sample_count = 0
    if cos_day_cycles != 1:
        while max(cos_cycle_index) >= cos_cycle_step:
            for spots in cos_cycle_index:
                if spots == cos_cycle_step:
                    sample_count += 1
                elif spots > cos_cycle_step:
                    cos_sample_time_per_day_list.append(sample_count)
                    sample_count = 1
                    cos_cycle_step += 1
                else:
                    cos_cycle_step += 1

    if cos_day_cycles != 1:
        cos_sample_time_per_day_list.append(sample_count)

    # generate dictionaries of lists each containing time & day cycles
    cos_daily_time = {}
    cos_daily_data = {}
    for day in range(0, cos_day_cycles):
        start = len(x) - sum(cos_sample_time_per_day_list[day:])
        end = len(x) - sum(cos_sample_time_per_day_list[day + 1:])
        cos_daily_time[day + 1] = [x[i] for i in range(start, end)]
        cos_daily_data[day + 1] = [y[i] for i in range(start, end)]

    # separate the real cycles from pseudo cycles, e.g. double dips
    real_onset_coords, real_offset_coords = get_real_cycles(cos_all_time, cos_onset_coords, cos_offset_coords, cos_fitted_on, cos_fitted_off, cos_day_cycles, cos_params[3])
    real_onsets = [cos_all_time[real_onset_coord] for real_onset_coord in real_onset_coords]
    real_offsets = [cos_all_time[real_offset_coord] for real_offset_coord in real_offset_coords]

    pseudo_onset_coords = [cos_onset_coord for cos_onset_coord in cos_onset_coords if cos_onset_coord not in real_onset_coords]
    pseudo_offset_coords = [cos_offset_coord for cos_offset_coord in cos_offset_coords if cos_offset_coord not in real_offset_coords]
    pseudo_onsets = [cos_all_time[pseudo_onset_coord] for pseudo_onset_coord in pseudo_onset_coords]
    pseudo_offsets = [cos_all_time[pseudo_offset_coord] for pseudo_offset_coord in pseudo_offset_coords]

    real_index_coords = sorted(real_onset_coords + real_offset_coords)
    cos_peak_time_list, cos_peak_data_list = peak_time_data(cos_all_time, cos_all_data, real_index_coords, real_onset_coords)

    # trapezoidal area under curve function to find auc between all qualifying onset and offset times
    # uses cos_peak_time_list, and cos_peak_data_list
    def peak_auc(time, data):
        if len(time) < 1:
            first_last_duration_amount = 0
            log_error("Peak not found.")
        elif (2*len(time))*(data[0]+data[len(data)-1]) != 0:
            first_last_duration_amount = ((time[len(time)-1] - time[0])/((2*(len(time)-1))*(data[0]+data[len(data)-1])))
        else:
            first_last_duration_amount = 0

        # all other indices
        tote = 0
        for pts in data:
            if pts != data[0] or pts != data[len(data)-1]:
                mid_duration_amount = ((time[len(time)-1]) - time[0])/(2*(len(time)-1))*2*pts
                tote += mid_duration_amount
        night_total = tote + first_last_duration_amount
        return night_total

    # write all trapezoidal auc computations to list called cos_peak_auc_list
    peak_cycle = 0
    cos_peak_auc_list = []

    while peak_cycle < len(cos_peak_time_list):
        cos_peak_auc_list.append(peak_auc(cos_peak_time_list[peak_cycle], cos_peak_data_list[peak_cycle]))
        peak_cycle += 1

    # write all peak_mp_auc calculations to a list
    mp_cycle = 0
    cos_peak_mp_auc_list = []

    while mp_cycle < len(cos_peak_time_list):
        cos_peak_mp_auc_list.append(midpoint_peak_auc(cos_peak_time_list[mp_cycle], cos_peak_data_list[mp_cycle]))
        mp_cycle += 1

    cos_peak_mp_auc_list_unstring = ", ".join("{0:.3f}".format(num) for num in cos_peak_mp_auc_list)
    cos_on_off_auc = (cos_peak_mp_auc_list[0] if len(cos_peak_mp_auc_list) == 1 else cos_peak_mp_auc_list_unstring)

    # setup the figure plots to be generated
    rcParams["figure.figsize"] = (10, 15)
    rcParams["legend.fontsize"] = 16
    rcParams["axes.labelsize"] = 16

    # create a blank figure
    fig = plt.figure()

    # add a plot (2x1 grid in the 1st position)
    axes = fig.add_subplot(211)

    # generate smooth fitted curves by upping the resolution to 100
    cos_time_fit = np.linspace(np.min(time), np.max(time), 100)
    cos_data_fit = cos_fit(cos_params, cos_time_fit)
    cos_acro_x_y = (cos_acro, cos_peak_value)
    # plot the data ("ro" = red circles) and the fit ("r-" = red line)
    axes.plot(time, data, "k-", label="Actual Data")
    axes.plot(time, data, "ko")
    axes.plot(cos_time_fit, cos_data_fit, "r-", label="Cosine Fit")
    axes.plot(cos_acro_list_x, cos_acro_list_y, "ro")
    plt.hlines(cos_mesor, time[0], time[len(time) - 1], "y", label="Threshold")
    # axes.plot(cos_crossing_points, cos_y_int, "yo")
    if cos_y_int:
        axes.plot(pseudo_onsets, [cos_y_int[0]] * len(pseudo_onsets), "yo")
        axes.plot(pseudo_offsets, [cos_y_int[0]] * len(pseudo_offsets), "yo")
        axes.plot(real_onsets, [cos_y_int[0]] * len(real_onsets), "o", color="orange")
        axes.plot(real_offsets, [cos_y_int[0]] * len(real_offsets), "o", color="orange")

    # add a legend
    axes.set_title("Ski Slope Cosine Fit" + " (" + name_of_dataset + ")")
    axes.set_xlabel("Hour")
    axes.set_ylabel("Measurement")
    axes.legend()

    # add a text box below the graph
    axes2 = fig.add_subplot(212)
    axes2.axis("off")

    # append additional calculations below the graph
    text = f"SS = {ss:,.3f}   RSD = {RSD:,.3f}\n\n"
    text += f"$r({res.size - 2:d}) = {r:,.3f}$,  $r^2({res.size - 2:d}) = {r2:,.3f}$\n\n"
    text += f"Peak Coordinates = ({cos_acro:,.3f}, {cos_peak_value:,.3f})\n\n"
    text += f"{ask_user_threshold:,.0f}% Threshold = {user_threshold:,.3f}   Number of Threshold Crossings = {len(cos_crossing_points):d}\n\n"
    text += f"Times of Threshold Crossings:\n{cos_crossing_points_unstring:s}\n\n"
    text += f"On-Off AUC:\n{cos_peak_mp_auc_list_unstring:s}\n\n"
    text += f"Fitted Onset = {cos_fitted_on:.3f}\nFitted Offset = {cos_fitted_off:.3f}\nFitted Duration = {cos_fitted_duration:.3f}\n\n"
    text += f"h = {cos_params[0]:.3f}, b = {cos_params[1]:.3f}, v = {cos_params[2]:.3f}, p = {cos_params[3]:.3f}"

    axes2.text(.02, .96, text, color="k", fontsize=16, horizontalalignment="left", verticalalignment="top", wrap=True,
               transform=axes2.transAxes)

    # save the figure as an image and close pyplot
    fig.savefig(destpath)
    plt.close()


if __name__ == "__main__":
    main()

