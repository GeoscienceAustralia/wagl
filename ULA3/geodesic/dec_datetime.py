#!/usr/bin/env python
'''
Copied across as-is from old ULA code.
'''
#TODO: Clean up code or replace
import datetime


def hour_fraction(dt):
    """
    Convert the minutes, seconds and microseconds of a :py:class:`datetime.datetime` object to a
    fraction of an hour.

    :param dt:
        Time to convert to a fraction of an hour.
    :type dt:
        :py:class:`datetime.datetime`

    :return:
        The minutes, seconds and microseconds converted to a fraction of an hour.
    """

    return (
        float(dt.minute) / 60 +
        float(dt.second) / 3600 +
        float(dt.microsecond) / 3600000000
    )


def day_fraction(dt):
    """
    Convert the hour, minutes, seconds and microseconds of a :py:class:`datetime.datetime` object to a
    fraction of a day.

    :param dt:
        Time to convert to a fraction of an hour.
    :type dt:
        :py:class:`datetime.datetime`

    :return:
        The hours, minutes, seconds and microseconds converted to a fraction of a day.
    """

    return (float(dt.hour) + hour_fraction(dt)) / 24


def decimal_day_of_month(dt):
    """
    Convert the day, hours, minutes, seconds and microseconds of a :py:class:`datetime.datetime` object to a
    decimal day of month.

    This is done by adding the day of the month to the fractional part of the day calculated from the
    hours, minutes, seconds and microseconds from ``dt``.

    :param dt:
        Time to convert to a decimial day.
    :type dt:
        :py:class:`datetime.datetime`

    :return:
        The days, hours, minutes, seconds and microseconds converted to a fraction of a day.
    """
    return float(dt.day) + day_fraction(dt)


def tu(dt):
    """
    Convert the hours, minutes, seconds and microseconds of a :py:class:`datetime.datetime` object to a
    decimal hour month.

    This is done by adding the day of the month to the fractional part of the day calculated from the
    hours, minutes, seconds and microseconds from ``dt``.

    :param dt:
        Time to convert to a decimial day.
    :type dt:
        :py:class:`datetime.datetime`

    :return:
        The days, hours, minutes, seconds and microseconds converted to a fraction of a day.
    """

    return float(dt.hour) + hour_fraction(dt)

