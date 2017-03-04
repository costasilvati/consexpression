# -*- coding: utf-8 -*-


class Message(object):
    """
    Contain all messages of system
    """

    def __init__(self):
        pass

    def message_1(self, append):
        """
        Valid name
        :param append:
        :return: void
        """
        print (" -- Valid " + append + "!")

    def message_2(self, append):
        """
        Error in experiment configuration
        :param append:
        :return: void
        """
        print ("********************************************\nError! The experiment don't have " + append + ".")

    def message_3(self, append):
        """
        uideline about corrections in configuration file
        :param append:
        :return: void
        """
        print ("=============================================\nCorrect the " +
               append + " in CONFIG_tool, and execute again.")
        print ("Close Session!")

    def message_4(self, append):
        """
        Default error message
        :param append:
        :return:
        """
        print ("Error! " + append + ".")

    def message_5(self, append):
        """
        File not found
        :param append:
        :return: void
        """
        print ("--> File not found: " + append + ".")

    def message_6(self, append):
        """
        Extension error
        :param append:
        :return: VOID
        """
        print ("Extension could be " + append + ".")

    def message_7(self, append):
        """
        Warning message
        :param append:
        :return: VOID
        """
        print ("Warning!! " + append + ".")

    def message_8(self, append):
        """
        Message whit header
        :param append:
        :return: VOID
        """
        print ('==============================\nMessage:\n ' + append + '.\n==============================')

    def message_9(self, append):
        """
        Only a message...
        :param append:
        :return:
        """
        print (append)