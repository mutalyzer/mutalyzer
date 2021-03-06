"""
Module for storing output and messages.

Output is stored as a named list that can be expanded. Messages can be
retrieved at a later time to provide flexibility. Message levels are
defined to increase or decrease the amount of logging and ouput.

The position of the log file, as well as the levels are defined in the
configuration file.

Message levels:
  - -1 : Log     ; Specifically log a message.
  -  0 : Debug   ; Debug information.
  -  1 : Info    ; Info.
  -  2 : Warning ; Regular warnings.
  -  3 : Error   ; Serious errors that can be compensated for.
  -  4 : Fatal   ; Errors that are not recoverable.
  -  5 : Off     ; Can be used as a log/output level to turn off output.

Public classes:
  - Message ; Container class for message variables.
  - Output  ; Output interface for errors, warnings and logging.
"""


from __future__ import unicode_literals

import io
import time

from mutalyzer import util
from mutalyzer.config import settings


class Output() :
    """
    Provide an output interface for errors, warnings and logging purposes.

    Private variables:
        - _outputdata ; The output dictionary.
        - _messages   ; The messages list.
        - _instance   ; The name of the module that made this object.
        - _loghandle  ; The handle of the log file.
        - _errors     ; The number of errors that have been processed.
        - _warnings   ; The number of warnings that have been processed.

    Special methods:
        - __init__(instance) ; Initialise the class with the calling
                               module.

    Public methods:
        - addMessage(filename, level, code, description) ; Add a message to
                                                           the message list.
        - getMessages()           ; Print all messages that exceed the
                                    configured output level.
        - addOutput(name, data)   ; Add output to the output dictionary.
        - getOutput(name)         ; Retrieve data from the output dictionary.
        - Summary()               ; Print a summary of the number of errors
                                    and warnings.
    """

    def __init__(self, instance):
        """
        Initialise the class with the calling module.

        Private variables (altered):
            - _outputdata ; The output dictionary.
            - _messages   ; The messages list.
            - _instance   ; Initialised with the name of the module that
                             created this object.
            - _loghandle  ; Initialised as the handle of the log file
                             defined in the configuration file.
            - _errors     ; Initialised to 0.
            - _warnings   ; Initialised to 0.

        @arg instance: The filename of the module that created this object
        @type instance: unicode
        """
        self._outputData = {}
        self._messages = []
        self._instance = util.nice_filename(instance)
        self._loghandle = io.open(settings.LOG_FILE, mode='a+',
                                  encoding='utf-8')
        self._errors = 0
        self._warnings = 0
    #__init__

    def addMessage(self, filename, level, code, description) :
        """
        Add a message to the message list.
        If the level exceeds the configured loglevel or if the level is -1,
        then the message is also logged.
        If the severity equals 2, then the number of warnings is inreased,
        if it exceeds 2, then the number of errors is increased.

        Private variables:
            - _messages  ; The messages list.
            - _instance  ; Module that created the Output object.
            - _loghandle ; Handle to the log file.

        Private variables (altered):
            - _warnings ; Increased by one if the severity equals 2.
            - _errors   ; Increased by one if the severity exceeds 2.

        @arg filename:    Name of the calling module
        @arg level:       Severity of the message
        @arg code:        Error code of the message
        @arg description: Description of the message
        """
        nice_name = util.nice_filename(filename)
        message = Message(nice_name, level, code, description)

        # Append a new message object to the messages list.
        self._messages.append(message)

        if level == 2:
            self._warnings += 1
        if level > 2:
            self._errors += 1

        # Log the message if the message is important enough, or if it is only
        # meant to be logged (level -1).
        if level >= settings.LOG_LEVEL or level == -1 :
            self._loghandle.write(time.strftime(
                settings.LOG_TIME_FORMAT + ' ') + "%s (%s) %s: %s: %s\n" % (
                    self._instance, nice_name, code, message.named_level(),
                    description))
            self._loghandle.flush()
        #if
    #addMessage

    def getMessages(self) :
        """
        Print all messages that exceed the configured output level.

        Private variables:
            - _messages  ; The messages list.

        @return: A list of messages
        @rtype: list
        """
        return filter(lambda m: m.level >= settings.OUTPUT_LEVEL,
                      self._messages)
    #getMessages

    def getMessagesWithErrorCode(self, errorcode):
        """
        Retrieve all messages that have a specific error code.

        Private variables:
            - _messages   ; The messages list.

        @arg errorcode: The error code to filter on
        @type errorcode: unicode

        @return: A filtered list
        @rtype: list
        """
        return filter(lambda m: m.code == errorcode, self._messages)
    #getMessagesWithErrorCode

    def getBatchMessages(self, level):
        """
        Returns a list of Messages with an errorlevel >= level
        and removes additional lines from a parseerror

        Private variables:
            - _messages   ; The messages list.

        @arg level: error level
        @type level: integer

        @return: list of Messages
        @rtype: list
        """
        ret = []
        lastorigin = ""
        for i in self._messages:
            if i.level >= level:
                # Todo: We changed this from 'Parser' to 'grammar', does this
                # still work?
                if lastorigin == 'grammar': #Only one parse error
                    continue
                lastorigin = i.origin
                ret.append("(%s): %s" % (i.origin, i.description))
            #if
        #for
        return ret
    #getBatchMessages

    def addOutput(self, name, data) :
        """
        If the output dictionary already has a node with the specified
        name, the list that this name points to is expanded with the data.
        Otherwise create a node and assign a list containing the data.

        Private variables:
            - _outputData ; The output dictionary.

        @arg name: Name of a node in the output dictionary
        @type name: unicode
        @arg data: The data to be stored at this node
        @type data: object
        """
        if self._outputData.has_key(name) :
            self._outputData[name].append(data)
        else :
            self._outputData[name] = [data]
    #addOutput

    def getOutput(self, name) :
        """
        Return a list of data from the output dictionary.

        Private variables:
            - _outputData ; The output dictionary.

        @arg name: Name of a node in the output dictionary
        @type name: string

        @return: output dictionary
        @rtype: dictionary
        """
        if self._outputData.has_key(name) :
            return self._outputData[name]
        return []
    #getOutput

    def getIndexedOutput(self, name, index, default=None):
        """
        Return an element of a list, the list is called 'name' in de
        _outputData dictionary. If either the list or the element does not
        exist, return {default}.

        @arg name:  Name of the list.
        @arg index: Index of the element to be retuned.
        @arg default: Default to return if either the list or the element
            does not exist.

        Private variables:
            - _outputData ; The output dictionary.

        @return: The requested element or None
        @rtype: any type
        """
        if self._outputData.has_key(name) :
            if 0 <= index < len(self._outputData[name]) :
                return self._outputData[name][index]
        return default
    #getIndexedOutput

    def Summary(self) :
        """
        Print a summary of the number of errors and warnings.

        Private variables:
            - _errors   ; The number of errors.
            - _warnings ; The number of warnings.

        @return:
            triple:
                - Number of errors
                - Number of warnings
                - Summary
        @rtype: integer, integer, unicode
        """
        e_s = 's'
        w_s = 's'
        if self._errors == 1 :
            e_s = ''
        if self._warnings == 1 :
            w_s = ''

        return self._errors, self._warnings, "%i Error%s, %i Warning%s." % (
            self._errors, e_s, self._warnings, w_s)
    #Summary
#Output

class Message() :
    """
    Container class for message variables.

    Special methods:
        - __init__(origin, level, code, description) ; Make a message object.

    Public variables:
        - origin      ; Name of the module creating this object.
        - level       ; Importance of the message.
        - code        ; The error code of the message.
        - description ; A description of the message.
    """

    def __init__(self, origin, level, code, description) :
        """
        Make a new message object.

        Public variables (altered):
            - origin      ; Name of the module creating this object.
            - level       ; Importance of the message.
            - code        ; The error code of the message.
            - description ; A description of the message.

        @arg origin: Name of the module creating this object
        @type origin: unicode
        @arg level: Importance of the message
        @type level: integer
        @arg code: The error code of the message
        @type code: unicode
        @arg description: A description of the message
        @type description: unicode
        """
        self.origin = origin
        self.level = level
        self.code = code
        self.description = description
    #__init__

    def __repr__(self):
        return 'Message("%s", %i, "%s", "%s")' % \
               (self.origin, self.level, self.code, self.description)
    #__repr__

    def __unicode__(self):
        return '%s (%s): %s' % \
               (self.named_level(), self.origin, self.description)
    #__unicode__

    def named_level(self):
        """
        Get message log level as readable string.

        @return:     A readable description of the log level.
        @rtype:      unicode
        """
        if self.level == 0:
            return "Debug"
        if self.level == 1:
            return "Info"
        if self.level == 2:
            return "Warning"
        if self.level == 3:
            return "Error"
        if self.level == 4:
            return "Fatal"
        return ''
    #named_level
#Message
