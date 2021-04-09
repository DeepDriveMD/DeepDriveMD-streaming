import datetime

def get_now():
    return datetime.datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
