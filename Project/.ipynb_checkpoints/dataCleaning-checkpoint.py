import pandas as pd
import numpy as np

"""
Data only record call option
Data should took the form of a dictionary
with keys = dates of buying
value = DataFrame that carries all info
"""


def append_risk_return(risk_data_name, placeholder):
    """
    :param risk_data_name: name of LIBOR U.S. dollar risk free interest rate data
    :param placeholder: dic:
    :return: dic: key are datetime object; values are risk
    """
    risk_data = pd.read_csv("./{}.csv".format(risk_data_name), index_col=0)
    risk_data.index = pd.to_datetime(risk_data.index)

    keys = list(placeholder.keys())
    temp = pd.DataFrame(data=np.array(keys), columns=["fuck"])
    ref = sorted(pd.to_datetime(temp["fuck"]))

    last_risk = 2
    last_price = 2680
    for i, p in enumerate(ref):
        if p in list(risk_data.index):
            if risk_data[risk_data_name].loc[p] != ".":
                last_risk = float(risk_data[risk_data_name].loc[p])/100
        placeholder[keys[i]][1] = last_risk
        curr_price = placeholder[keys[i]][0]["spotclose"].values[0]
        placeholder[keys[i]][2] = (curr_price - last_price) / last_price
        last_price = curr_price


class DataCleaner:
    def __init__(self, columns, date_columns, verbose=True, call_only=True):
        """
        :param columns: a list of data required fields
        :param date_columns: a list of date related
        """
        self.verbose = verbose
        self.columns = columns
        self.dateColumns = date_columns
        self.num_days = 0
        self.num_observation = 0
        self.callOnly = call_only
        if self.verbose:
            print("Required fields are:{}".format(str(columns)))

    def trim(self, data_name, top):
        """
        :param data_name: a string of datafile name
        :param top: the # of mostly traded ticker on the day
        :return: tuple, (date, trimmed pandas dataFrame)
        """
        month_data = pd.read_csv("./SPXdata/{}.csv".format(data_name))[self.columns]

        # Check if we require only call option
        if self.callOnly:
            month_data = month_data[month_data["OptionType"] == "call"]

        # transform date columns into date type
        for p in self.dateColumns:
            month_data[p] = pd.to_datetime(month_data[p].astype(str), format="%Y%m%d")

        # create time to exercise column
        month_data["T"] = month_data["expirydate"] - month_data["date"]

        # a list of trading days datetime objects in the every month
        keys = month_data["date"].unique()

        # a list of trading days data in the dataset
        values = []

        day_delta = 0
        observation_delta = 0
        for i, p in enumerate(keys):
            day_data = month_data[month_data["date"] == p]
            values.append(day_data.nlargest(top, "volume"))
            day_delta += 1
            observation_delta += top
        del month_data

        # verbose
        if self.verbose:
            print("processed {} days".format(day_delta))
            print("processed {} observations".format(observation_delta))

        # add processed data info
        self.num_days += day_delta
        self.num_observation += observation_delta

        return keys, values

    def insert_data(self, placeholder, data_name, risk_data_name, top):
        """
        :param placeholder: dictionary that holds the final data
        :param data_name: a string of datafile name, no .csv at last
        :param risk_data_name: name of LIBOR U.S. dollar risk free interest rate data
        :param top: # of mostly traded ticker
        :return: void
        """
        days, records = self.trim(data_name, top)
        for i in range(len(days)):
            placeholder[days[i]] = [records[i], None, None]
        append_risk_return(risk_data_name, placeholder)

    def get_num_days(self):
        print("Attached {} days of data".format(self.num_days))
        return self.num_days

    def get_num_observation(self):
        print("Attached {} option observations".format(self.num_observation))
        return self.num_observation


# if __name__ == '__main__':
#     cleaner = DataCleaner(columns=["code", "OptionType", "expirydate", "date", "strike", "spotclose"],
#                           date_columns=["date", "expirydate"])
#     data = {}
#     cleaner.insert_data(data, "SPX_201704_1day", 3)
#
#     cleaner.append_risk("USD3MTD156N", data, "4/3/2017", "4/28/2017")
