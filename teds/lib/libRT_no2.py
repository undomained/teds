import datetime
import numpy as np
import time
import subprocess
import sys
import logging
import inspect
import os
import sys
import json
import tables
import warnings
import re
from collections import OrderedDict
import tempfile
import shutil


##### disamar tools


# rt_run and rt_cfg from: 

# Maarten Sneep, 2009
# Purpose: parser for a configuration file for the radiative transfer module
#

'''
DISAMAR Configuration file parser
=================================

  Read and write a configuration file for Johan's DISIMAR program.

  The format Johan suggested is somewhat convoluted, with sections,
  subsections, and the possibility of repeated keys within a subsection.
  The repeated keys give values for subsequent channels, so order is
  very important. To add to that, a key can contain multiple values
  on a single line. All parts are separated by whitespace.

  Comments are indicated by a # mark at the start of a line, or between ()
  after the last element on a line, and should be preserved as best we can.

  Using '&' at the end of a line, indicates that the next line continues
  the current one.

  In normal operation repeated section, subsections and keys add to,
  rather than replace existing sections. This can be controlled by adding a
  meta-comment. With the following line
  #@ replace_section = True
  existing sections will be replaced from that point in the file onwards.
  The same can be limited to subsections::

    #@ replace_subsection = True

  And finally, keys can be replaced, rather than added to with::

    #@ replace_keys = True

  Note that in the latter mode, the original behavior is I{severely} altered.

  In short: the format is almost too flexible. An example is included at
  the end of this doc string.

  Note that Johan has extensive testing of the contents of the configuration
  file. This class does not attempt to duplicate that behaviour. As long as
  the format of configuration is technically correct, it will be read and
  declared valid.

  Limitations
  -----------

    - Single line # comments are moved to the top of the file on output
    - The order in which key are emitted is random, only the relative order
      of repeated keys is guaranteed.

  @author: Maarten Sneep
  @contact: maarten.sneep@knmi.nl
  @organization: KNMI
  @copyright: Maarten Sneep, KNMI
  @license: Creative Commons Attribution-ShareAlike 3.0 Netherlands License.

  Example file
  ------------

  A short section of a sample configuration file::

    SECTION GENERAL
    subsection fitAltAerCldSurf
    fit_cloud_top         1
    fit_cloud_fraction    0

    subsection filenames
    referenceSpecFileName  refSpecHighLat.dat   (radiance reference spectrum)
    polCorrectionFileName  polCorLat-45.dat     (polarization correction file)

    SECTION INSTRUMENT_SPECIFICATIONS
    subsection FWHM (in nm)
    FWHM_Irradiance_sim    0.45
    FWHM_radiance_sim      0.45d0
    FWHM_Irradiance_retr   0.45
    FWHM_radiance_retr     0.45
    FWHM_Irradiance_sim    0.35
    FWHM_radiance_sim      0.35e0
    FWHM_Irradiance_retr   0.35
    FWHM_radiance_retr     0.35

    subsection slit_index
    # single line comment
    slit_index_Irradiance_sim  1
    slit_index_radiance_sim    1
    slit_index_Irradiance_retr 1
    slit_index_radiance_retr   1

    SECTION SCENE_PARAMETERS
    subsection albedo
    surfacealbedo 0.05 0.04 0.03 0.02 0.02

  Full configuration files are included in the C{InputFiles} dicrecory in the
  DISAMAR distribution.

  Usage
  -----

  A short sample on how to use the code::

    # load the package
    import rt_cfg

    # create configuration object, specifying a file will cause it to be read.
    cfg = rt_cfg.RT_configuration(filename='test.config')

    # select a configuration item
    item = cfg['GENERAL','filenames','referenceSpecFileName']

    # item is now a configuration item object. Its value can be obtained with
    print(item.values()) # which prints 'refSpecHighLat.dat'

    # change the value of an item
    cfg['SCENE_PARAMETERS','albedo','surfacealbedo'].setvalue([0.05, 0.04, 0.04])

    # Add an element. If you want to add a comment, use a tuple.
    cfg['SCENE_PARAMETERS','albedo','cloud'] = (0.8, 'comment')

    # change the file name
    cfg.setfile('new.config')

    # write to the new file
    cfg.write()
'''

# import integer division with sane results


if sys.version_info[0] == 3:
    string_types = (str,)
else:
    string_types = (str,)


def dblformat(item):
    if type(item) == type(1.0):
        s = str(item).lower()
        if 'e' in s:
            rval = s.replace('e', 'D')
        else:
            rval = "{0}D0".format(s.strip())
    else:
        rval = str(item)
    return rval

class RT_dict(OrderedDict):
    '''
An extended dictionary

This extended dictionary includes a comment-field and a comment-block. The
base-class is an ordered dictionary, standard in python 2.7, but defined above
for older python versions. The contents are other L{RT_dict} objects or
L{RT_item} objects.

@ivar _comment: Comment for the section or subsection
@type _comment: C{string}
@ivar _blockcomment: Comment block for the file, section or subsection
@type _blockcomment: C{string}
'''
    def __init__(self, *args, **kwargs):
        '''Initializer

@param args: Positional arguments for the C{dict} constructor
@type args: C{[]}
@param kwargs: Keyword parameters for the C{dict} constructor. A few are
               filtered for use as comments.
@type kwargs: C{{}}

@keyword comment: A comment for the section or subsection, the one between parentheses
@keyword blockcomment: Comments for the block, inserted as #-marked lines.

The specified keyword are removed after they are used. The rest of the arguments
(positional and keywords) is passed on to C{super()}.

@note: This class is now much simpler because the OrderedDict is used from Python 2.7
'''
        if 'comment' in kwargs:
            comment = kwargs['comment']
            del kwargs['comment']
        else:
            comment = ''

        if 'blockcomment' in kwargs:
            blockcomment = kwargs['blockcomment']
            del kwargs['blockcomment']
        else:
            blockcomment = ''

        if len(args) >= 1:
            super(self.__class__, self).__init__(args[0])
        else:
            super(self.__class__, self).__init__(self)

        self._comment = comment
        self._blockcomment = blockcomment

    def setcomment(self, comment):
        '''Change the comment for the object

@param comment: The new comment for the group.
@type comment: C{string}.
'''
        self._comment = comment

    def comment(self):
        '''Obtain the current comment

@return: The current comment.
@rtype: C{string}
'''
        return self._comment

    def setblockcomment(self, comment):
        '''Build a longer comment block

The comments are appended to any existing block comments, separating them with a newline and a hash-mark.

@param comment: The new comment for the group.
@type comment: C{string}.
'''
        if not self._blockcomment:
            self._blockcomment = comment
        else:
            self._blockcomment = self._blockcomment + '\n#' + comment

    def blockcomment(self):
        '''Obtain the current block comment

@return: The current block comment,.
@rtype: C{string} (newline and hash-mark separated).
'''
        return self._blockcomment

    def add_item(self, item):
        self[item.key()] = item

    def __str__(self):
        return self.blockcomment() + "\n" + "\n".join([str(v) for v in list(self.values())])

class RT_item(object):
    '''A configuration item.

This object contains the name of the item, the value(s) and the comment string.
This is the finest granularity of a configuration file.

@ivar _value: The container for the values, most likely a list, or even a list of lists.
@ivar _comment: The parenthesized comment string for the item.
@ivar _key: The key to identify the object.
'''

    def __init__(self, key=None, value=None, comment=None):
        '''Initializer

@param value:   can be an existing 'item' object, in which case the remaining
         parameters are ignored. If it is a string, it will be parsed,
         otherwise the value is taken as-is.
@param key:     Required if value is not an 'item' object. String.
@param comment: A comment string for the item.
'''
        self._value = []
        self._comment = []
        if type(value) == type(self):
            self._key = value.key()
            self.setvalue(list(value.values()))
            self._comment = value.comments()
        else:
            if key is None:
                raise ValueError('a key is required')
            self._key = key
            if type(value) == type(''):
                self.appendvalue(value)
            else:
                self.setvalue(value)
            self.appendcomment(comment)

    def convertvalue(self, value): # IGNORE:R0201
        '''Convert a value as read from file to int, float, or string [list].

Always returns a list of item(s), even if there is a single value.
This is convenient when writing the file later. The converter splits the input
string on whitespace. For each item it first tries  to convert the string to an
integer. If that fails, it tries to convert the value to a float value. If that
fails it uses a string value as a fallback option.

@param value: The raw input as read from the file.
@type value: C{string}
'''
        # test if we have a value at all
        if value is None:
            value = ""

        # Split on whitespace
        if isinstance(value, str):
            vallist = value.split()
        elif isinstance(value, list):
            vallist = [str(item) for item in value]
        elif isinstance(value, (float, int)):
            vallist = [str(value)]
        else:
            raise TypeError('Unexpected type for value (got "{0}")'.format(type(value)))

        newvallist = []

        # loop over the elements of the value string
        for elem in vallist:
            # convenience block to try each conversion.
            # Start with the most stringent type (int),
            # and progress towards the ost relaxed (string, no conversion needed)
            while True:
                try:
                    newval = int(elem)
                    break
                except ValueError:
                    pass

                try:
                    # the float interpreter is C-based, and does not like 1d0 as an input
                    newval = float(elem.lower().replace('d', 'e'))
                    break
                except ValueError:
                    pass

                newval = elem
                break
            # Add the converted value to the output list
            newvallist.append(newval)
        return newvallist

    def append(self, key, value, comment):
        '''Append a value and a comment to the item

@param key: The key is included in the call for sanity checking. If the provided
            key is not equal to our stored key, then we fail.
@param value: The new value. This value is added with the L{appendvalue} method.
@param comment: The new or additional comment for the key/value pair. This
            comment is added with the L{appendcomment} method.
'''
        # check that we refer to the same key
        if key != self._key:
            raise ValueError(format('The given key ({0}) does not match the stored key ({1})',
                             key, self.key))
        self.appendvalue(value)
        self.appendcomment(comment)

    def appendvalue(self, value):
        '''Append a value to the item.

The value is first converted with the L{convertvalue} method.

@param value: The raw (string) input.
'''
        newval = self.convertvalue(value)
        if self._value is None:
            self._value = []
        self._value.append(newval)

    def appendcomment(self, comment):
        '''Append a comment to the list of comments.

A value is always appended to the list of comments, even if it is just the empty string.

@param comment: The new comment (can be C{None}).
'''
        if comment is None or comment == 'None':
            comment = ''
        if self._comment is None:
            self._comment = []
        if type(comment) == type(""):
            comment = [comment]

        if type(comment) == type([]):
            self._comment = comment[:]
            new_comment_length = len(self._value)
            old_comment_length = len(self._comment)
            if old_comment_length < new_comment_length:
                for dummy in range(new_comment_length - old_comment_length):
                    self._comment.append('')
            else:
                self._comment = self._comment[0:new_comment_length]

    def comments(self):
        '''Retrieve the comments for this item.

@return: The list of comments (shallow copy of the internal list).
@rtype: C{[]}
'''
        if len(self._comment) > 1:
            return self._comment[:]
        else:
            return self._comment[0]

    def values(self):
        '''Retrieve the values for this item.

The values are always stored as a list of lists. Here we unpack the
items as far as possible.

@return: The value of the item.
@rtype: atom, list, or list of lists.
'''
        if len(self._value) == 1:
            if len(self._value[0]) == 1:
                return self._value[0][0]
            else:
                return self._value[0]
        else:
            return self._value[:]

    value = values

    def rawvalues(self):
        '''get the raw value. Makes it easier to loop over multiple bands.

Returns a shallow copy of the values list.

@return: values.
@rtype: C{[[]]}
'''
        return self._value[:]

    def set_rawvalue(self, newvalue):
        '''Set a raw value.

Ensure that the comments list has the same length as the values list. Otherwise
hardly any checks are performed on the input parameter. Caveat emptor!

@param newvalue: The new (raw) value to use.
@type newvalue: C{[[]]}
'''
        self._value = newvalue
        if len(self._value) != len(self._comment):
            new_comment_length = len(self._value)
            old_comment_length = len(self._comment) if self._comment is not None else 0
            if old_comment_length < new_comment_length:
                for dummy in range(new_comment_length - old_comment_length):
                    self.appendcomment('')
            else:
                self._comment = self._comment[0:new_comment_length]

    def setvalue(self, newvalue):
        '''Change the value for this key.

Ensure that the comments list has the same length as the values list. The values
are checked to ensure that they are packaged eventually as a list of lists of atoms.
A new 'raw' value is produced, and L{set_rawvalue} is used to update the value of C{self}.

@param newvalue: The new (raw) value to use.
@type newvalue: C{[[]]}
'''

        if type(newvalue) != type([]):
            rawvalue = [[newvalue]]
        else:
            if type(newvalue[0]) != type([]):
                rawvalue = [newvalue]
            else:
                rawvalue = newvalue
        self.set_rawvalue(rawvalue)


    def key(self):
        '''Retrieve the key for this item

@return: key
@rtype: C{string}
'''
        return self._key

    def extend(self, valstr, comment):
        '''Continue a configuration line, and extend the last item.

The format of the configuration file allows for continued lines. In that case the
last item of the values list must be extended, rather than a new list appended.

@param valstr: A raw string value.
@type valstr: C{string}
@param comment: Comment for this value
@type comment: C{string}
'''
        values = self.convertvalue(valstr)
        self._value[-1].extend(values)
        if comment:
            if self._comment[-1]:
                self._comment[-1] = self._comment[-1] + '; ' + comment
            else:
                self._comment[-1] = comment

    def __repr__(self):
        '''Create a representation of the values of this item.

The representation is suitable recreating self in Python.
'''
        vals = self.values()
        try:
            _ = (e for e in vals)
            if isinstance(vals, str):
                vals = [vals]
            else:
                vals = list(vals)
        except TypeError:
            vals = [vals]
        return '{cls}(key={key!r}, \n\tvalue={val!r}, \n\tcomment={comment!r})'.format(
                cls=self.__class__.__name__,
                key=self.key(),
                val=vals,
                comment=self.comments())

    def __str__(self):
        '''Produce a string representation in the correct syntax for the configuration file format.

The string output is suitable for generating an item in a configuration file.
'''
        out = []
        for val, comment in zip(self._value, self._comment):
            if comment:
                out.append('{0}  {1}  ({2})'.format(
                           self._key, ' '.join([dblformat(item) for
                           item in val]), comment))
            else:
                out.append('{0}  {1}'.format(self._key, ' '.join([dblformat(item) for
                                             item in val])))
        return '\n'.join(out)


class RT_configuration(object):
    '''Radiative transfer configuration manager and file parser

This is the main entrypoint for this module.'''

    def __init__(self, filename=None, config=None, comment=None, debug=False):
        '''Initializer

Parameters:
        filename:    optional, configuration filename.
        config:  optional, existing RTParser object.
        comment: optional, a comment at the filename-level.'''
        self._file = filename
        if config is not None and type(config) == type(self):
            self.config = config.dictionary()
        else:
            self.config = RT_dict(comment=comment)

        self.logger = logger()
        if debug:
            self.logger.setLevel(10)
        self.logger.debug("Initializing {0}".format(self.__class__.__name__))

        if not (self._file is None or (config is not None and type(config) == type(self))):
            self.logger.info("Reading config from {0}.".format(self._file))
            self.read()
        self.logger.debug("Created {0}".format(self.__class__.__name__))

    def dictionary(self):
        '''Return a copy of the underlying dictionary'''
        return self.config.copy()

    def file(self):
        '''Retrieve the file that the object is currently linked to.'''
        return self._file

    def setfile(self, filename):
        '''Set the file that will be used for reading or writing.'''
        self._file = filename

    def sections(self):
        '''Obtain a list of sections'''
        sections = list(self.config.keys())
        self.logger.debug("Sections in object: '{0}'.".format("', '".join(sections)))
        return sections

    def subsections(self, section):
        '''Obtain a list of subsections'''
        subsections = list(self[section].keys())
        self.logger.debug("Subsections in section '{1}': '{0}'.".format("', '".join(subsections), section))
        return subsections

    def keys(self, section, subsection):
        '''Obtain a list of keys in a particular subsection'''
        keys = list(self[section, subsection].keys())
        self.logger.debug("Keys in section '{1}', subsection {2}: '{0}'.".format("', '".join(keys), section, subsection))
        return keys

    def items(self, section, subsection):
        '''Obtain a list of key-value pairs for a particular subsection'''
        items = list(self[section, subsection].items())
        self.logger.debug("Items in section '{1}', subsection {2}: {{{0}}}.".format(
                ", ".join(["'{0}': '{1}'".format(k,v) for k,v in items]), section, subsection))
        return items

    def values(self, section, subsection):
        '''Return the values in the subsection'''
        self.logger.debug("Values for SECTION '{0}', subsection '{1}'".format(section, subsection))
        val = []
        for key in list(self[section, subsection].keys()):
            val.append(self[section, subsection, key])
        return val

    def __contains__(self, item):
        rval = 0
        if type(item) == type(''):
            rval = 1 if item in self.config else 0
        elif len(item) == 1:
            rval = 1 if item[0] in self.config else 0
        elif len(item) == 2:
            rval = 1 if item[0] in self.config and item[1] in self.config[item[0]] else 0
        else:
            rval =  1 if ((item[0] in self.config) and
                (item[1] in self.config[item[0]]) and
                (item[2] in self.config[item[0]][item[1]])) else 0
        self.logger.debug('__contains__({0}): {1}'.format(item, 'True' if rval else 'False'))
        return rval

    def __delitem__(self, key):
        '''Remove sections, subsections, and keys'''
        self.logger.debug('__delitem__: {0}'.format(key))
        if type(key) == type(''):
            del self.config[key]
        elif len(key) == 1:
            del self.config[key[0]]
        elif len(key) == 2:
            del self.config[key[0]][key[1]]
        else:
            del self.config[key[0]][key[1]][key[2]]

    def __getitem__(self, key):
        '''Get an item from the configuration dictionary.

Key can be a string or a tuple.
If value is a tuple, then it is assumed that the real value is the first element,
and the second element contains a comment'''
        self.logger.debug('__getitem__: {0}'.format(key))
        if isinstance(key, string_types):
            return self.config[key]
        elif len(key) == 1:
            return self.config[key[0]]
        elif len(key) == 2:
            return self.config[key[0]][key[1]]
        else:
            return self.config[key[0]][key[1]][key[2]]

    def __setitem__(self, key, value):  #IGNORE:R0915 #IGNORE:R0912
        '''Insert an item into the configuration.

Key can be a string (insert a whole section), or a tuple of strings
(insert a subsection, or a key into a subsection). If a section is already
there, then new subsections are merged, same for subsections and keys.

If a key is already there, then the new items extend the existing list,
or the single item is turned into a list, and the new value is appended.'''
        self.logger.debug('__setitem__: [{0}] = {1}'.format(key, value))
        # split the value into a value and a comment
        if type(value) == type(tuple()):
            comment = str(value[1])
            value = value[0]
        else:
            comment = ''

        # handle string or single item tuple values as key.
        # The value must be a dictionary describing the whole section
        if type(key) == type(tuple()) and len(key) == 1:
            key = key[0]

        if type(key) == type(''):
            if type(value) == type({}) or type(value) == type(RT_dict()):
                self.logger.debug("Inserting whole section at '{0}'.".format(key))
                if key in self.config:
                    self.config[key].update(value)
                else:
                    if value:
                        self.config[key] = RT_dict(value, comment=comment)
                    else:
                        self.config[key] = RT_dict(comment=comment)
                return
            else:
                self.logger.error('Value must be a dictionary if used in this way')
                raise ValueError('Value must be a dictionary if used in this way')

        # Key is a tuple with two items. Insert a whole subsection at once.
        if len(key) == 2:
            self.logger.debug("Inserting whole section at '{0[0]}.{0[1]}'.".format(key))
            if type(value) == type({}) or type(value) == type(RT_dict()):
                if key[0] not in self.config:
                    self.config[key[0]] = RT_dict()
                if key[1] in self.config[key[0]]:
                    self.config[key[0]][key[1]].update(value)
                else:
                    if value:
                        self.config[key[0]][key[1]] = RT_dict(value, comment=comment)
                    else:
                        self.config[key[0]][key[1]] = RT_dict(comment=comment)
                return
            else:
                self.logger.error('Value must be a dictionary if used in this way')
                raise ValueError('Value must be a dictionary if used in this way')

        # Key is a tuple with three items. This fully qualifies a single key item
        if len(key) == 3 and (type(value) != type({})
                         and type(value) != type(RT_dict())):
            self.logger.debug("Inserting item at '{0[0]}.{0[1]}.{0[2]}'.".format(key))
            if key[0] not in self.config:
                self.logger.debug("Creating section '{0[0]}'.".format(key))
                self.config[key[0]] = RT_dict()
            if key[1] not in self.config[key[0]]:
                self.logger.debug("Creating subsection '{0[0]}.{0[1]}'.".format(key))
                self.config[key[0]][key[1]] = RT_dict()
            if key[2] in self.config[key[0]][key[1]]:
                self.logger.debug("Appending to value '{0[2]}'.".format(key))
                self.config[key[0]][key[1]][key[2]].append(key[2], value, comment)
            else:
                self.logger.debug("Creating item '{0[0]}.{0[1]}.{0[2]}'.".format(key))
                self.config[key[0]][key[1]][key[2]] = RT_item(key=key[2],
                                                              value=value, comment=comment)
        else:
            self.logger.error('Value can not be a dictionary if used in this way')
            raise ValueError('Value can not be a dictionary if used in this way')

    def __str__(self):
        '''Create a string representation for writing to file.'''
        out = []

        self.logger.debug("Entering {0}.{1}".format(self.__class__.__name__, whoami()))
        file_comment = self.config.blockcomment()
        if file_comment:
            out.append('# {0}\n'.format(file_comment))

        # obtain the sections
        sections = self.sections()
        self.logger.debug(", ".join(sections))

        # loop over the sections
        for section in sections:
            # write the section header
            sect = self[section]
            comment = sect.comment()
            if comment:
                out.append('\n\nSECTION {0} ({1})\n'.format(section, comment))
            else:
                out.append('\n\nSECTION {0}\n'.format(section))

            comment = sect.blockcomment()
            if comment:
                out.append('# {0}\n'.format(comment))

            # obtain the subsections and put them in alphabetical order
            subsections = self.subsections(section)

            # loop over the subsections
            for subsection in subsections:
                # write the subsection header
                subsect = self[section, subsection]
                comment = subsect.comment()
                if comment:
                    out.append('\nsubsection {0} ({1})\n'.format(subsection, comment))
                else:
                    out.append('\nsubsection {0}\n'.format(subsection))
                comment = subsect.blockcomment()
                if comment:
                    out.append('# {0}\n'.format(comment))
                # loop over the keys in this subsection
                for key in list(self[section, subsection].keys()):
                    out.append("{0!s}\n".format(self[section, subsection, key]))
        self.logger.debug("Exit {0}.{1}".format(self.__class__.__name__, whoami()))
        return ''.join(out)

    def write(self, filename=None):
        '''
        Write self to file.

        If filename is given, then the file attribute is not used.
        If the file attribute of the is not valid and filename is
        not given, write to stdout.
        '''
        if filename is not None:
            fp = open(filename, 'w')
        else:
            if not self._file:
                fp = sys.stdout
            else:
                # Open the file for writing
                fp = open(self._file, 'w')

        fp.write(str(self))
        fp.write('\n')

        if self._file is not None or filename is not None:
            fp.close()

    def read(self, filename=None,
             replace_section=False, #IGNORE:R0915 #IGNORE:R0912
             replace_subsection=False,
             replace_keys=False):
        '''Read configuration from file (either the class attribute or given explicitly) or stdin'''
        if hasattr(self._file, "read"):
            fp = self._file
            lines = fp.readlines()
        elif filename and os.access(filename, os.F_OK) and os.access(filename, os.R_OK):
            # Open the file for reading
            fp = open(filename, 'r')
            lines = fp.readlines()
        elif self._file and os.access(self._file, os.F_OK) and os.access(self._file, os.R_OK):
            fp = open(self._file, 'r')
            lines = fp.readlines()
        elif isinstance(filename, string_types):
            # assume not a file is given but a string
            lines = filename.split('\n')
        else:
            fp = sys.stdin
            lines = fp.readlines()

        # define and compile regular expressions
        re_control = re.compile(
            r'''^#@[ \t]*([a-zA-Z_][a-zA-Z0-9_]+)[ \t]*=[ \t]*(.+?)$''')
        re_section = re.compile(
            r'''^[ \t]*SECTION[ \t]+([a-zA-Z_][a-zA-Z0-9_-]+)[ \t]*\(?(.*?)\)?[ \t]*$''')
        re_subsection = re.compile(
            r'''^[ \t]*subsection[ \t]+([a-zA-Z_][a-zA-Z0-9_-]+)[ \t]*\(?(.*?)\)?[ \t]*$''')
        re_key = re.compile(
            r'''^[ \t]*([a-zA-Z_][a-zA-Z0-9_-]*)[ \t]+([^(&]+)(?:\((.*?)\))?[ \t]*(\&)?[ \t]*$''')  #IGNORE:C0301
        re_continue = re.compile(
            r'''^[ \t]*([^(&]+)(?:\((.*?)\))?[ \t]*(\&)?[ \t]*$''')
        re_comment = re.compile(
            r'''^[ \t]*# ?(.*)$''')
        re_empty = re.compile(
            r'''^[ \t\n\r]*$''')
        re_linebreaks = re.compile(r'''[\r\n]{1,2}''')

        # set some starting values
        section = ''
        subsection = ''
        key = ''
        next_continues = False
        linecount = 0

        # loop over all lines
        for line in lines:
            linecount += 1
            try:
                # strip away linebreaks
                line = re_linebreaks.sub('', line)
                search = re_control.search(line)
                if search:
                    self.logger.debug('Line: {0} ("{1}") matched re_control'.format(linecount, line.replace('\n', '').replace('\r', '')))
                    if search.group(1).lower() == 'replace_section':
                        replace_section = (search.group(2).lower() == 'true')
                    elif search.group(1).lower() == 'replace_subsection':
                        replace_subsection = (search.group(2).lower() == 'true')
                    elif search.group(1).lower() == 'replace_keys':
                        replace_keys = (search.group(2).lower() == 'true')

                    if section:
                        if subsection:
                            self[section, subsection].setblockcomment('@ ' +
                                 search.group(1) + ' = ' + search.group(2))
                        else:
                            self[section].setblockcomment('@ ' +
                                 search.group(1) + ' = ' + search.group(2))
                    else:
                        self.config.setblockcomment('@ ' +
                             search.group(1) + ' = ' + search.group(2))
                    continue
                search = re_empty.search(line)
                if search:
                    self.logger.debug('Line: {0} matched re_empty'.format(linecount))
                    continue
                search = re_section.search(line)
                if search:
                    self.logger.debug('Line: {0} ("{1}") matched re_section'.format(linecount, line))
                    section = search.group(1)
                    subsection = None
                    if replace_section and section in self:
                        del self[section]
                    self[section] = ({}, search.group(2))
                    continue
                search = re_subsection.search(line)
                if search:
                    self.logger.debug('Line: {0} ("{1}") matched re_subsection'.format(linecount, line))
                    subsection = search.group(1)
                    if replace_subsection and (section, subsection) in self:
                        del self[section, subsection]
                    self[section, subsection] = ({}, search.group(2))
                    continue
                search = re_comment.search(line)
                if search:
                    self.logger.debug('Line: {0} ("{1}") matched re_comment'.format(linecount, line))
                    if section:
                        if subsection:
                            self[section, subsection].setblockcomment(' ' + search.group(1))
                        else:
                            self[section].setblockcomment(' ' + search.group(1))
                    else:
                        self.config.setblockcomment(' ' + search.group(1))
                    continue
                search = re_key.search(line)
                if search:
                    self.logger.debug('Line: {0} ("{1}") matched re_key'.format(linecount, line))
                    key = search.group(1)
                    valstr = search.group(2)
                    comment = search.group(3)
                    next_continues = bool(search.group(4))
                    if replace_keys and (section, subsection, key) in self:
                        del self[section, subsection, key]
                    self[section, subsection, key] = (valstr, comment)
                    continue
                search = re_continue.search(line)
                if search and next_continues:
                    self.logger.debug('Line: {0} ("{1}")\nmatched re_continue'.format(linecount, line))
                    valstr = search.group(1)
                    comment = search.group(2)
                    next_continues = bool(search.group(3))
                    self[section, subsection, key].extend(valstr, comment)
                    continue
                raise ValueError('Encountered some unexpected syntax at line # {1}:\n\t{0}'.format(
                                 line, linecount))
            except:
                print(("configline ({0:d}):\n{1}".format(linecount, line)))
                raise

        if (self._file and os.access(self._file, os.F_OK) and os.access(self._file, os.R_OK) or
            filename and os.access(filename, os.F_OK) and os.access(filename, os.R_OK)):
            fp.close()








class rt_run(object):  #IGNORE:C0103
    '''run the radiative transfer model (and retrieval system) DISAMAR
    '''
    def __init__(self, cfg=None,
                 disamar=None,
                 output=None,
                 spectrum=None,
                 debug=False,
                 quiet=False,
                 tempbase=None):
        '''@type cfg: rt_cfg.RT_configuration object or path to config file.
@type disamar: string (path to DISAMAR executable)
@type output: filename for output.
'''
        if type(cfg) == type(""):
            if os.path.exists(cfg):
                self._cfgfile = os.path.abspath(cfg)
                self.config = RT_configuration(self._cfgfile)
            else:
                raise ValueError("File {0} not found".format(cfg))
        else:
            self._cfgfile = os.path.abspath(cfg.file())
            self.config = cfg

        self._debug = debug
        self._quiet = quiet

        if not disamar:
            disamar = "/usr/people/sneep/Documents/Disamar/Disamar.exe"
        if not os.path.exists(disamar):
            raise ValueError('Disamar executable not found: {0}'.format(disamar))
        self.disamar = disamar

        if spectrum is not None and os.path.exists(spectrum):
            self.spectrum = os.path.abspath(spectrum)
        else:
            self.spectrum = None

        if output:
            self.output = os.path.abspath(output)
        else:
            f = self._cfgfile
            self.output = os.path.join(os.path.dirname(f), os.path.splitext(os.path.basename(f))[0] + '.hdf5')

        self.attributes = None

        fmt = "{time[0]:04d}{time[1]:02d}{time[2]:02d}_{time[3]:02d}{time[4]:02d}{time[5]:02d}_{pid:05d}"
        date_string = fmt.format(time=time.gmtime()[0:6], pid=os.getpid())

        if tempbase is not None and os.path.exists(tempbase):
            tempfile.tempdir = os.path.abspath(tempbase)

        self.tmploc = os.path.join(tempfile.gettempdir(), 'DISAMAR.' + date_string)
        os.mkdir(self.tmploc)
        self.make_copy()

    def __del__(self):
        if hasattr(self, 'tmploc'):
            if hasattr(self, '_debug') and not self._debug:
                shutil.rmtree(self.tmploc, True)

    def save_config(self):
        '''Save the configuration object in file 'fname'. '''
        self.copy_spectrum()
        self.config.write(filename=os.path.join(self.tmploc, 'Config.in'))

    def copy_spectrum(self):
        if self.spectrum and os.path.exists(self.spectrum):
            if self._debug:
                print('using spectrum from {0}'.format(self.spectrum))
            filename = os.path.basename(self.spectrum)
            self.config['GENERAL','overall','useReflFromFile'].setvalue(1)
            self.config['GENERAL','fileNames','externalReflectanceFileName'].setvalue(filename)
            os.symlink(self.spectrum, os.path.join(self.tmploc, filename))
        else:
            if self._debug:
                print('Simulate spectrum')
            self.config['GENERAL','overall','useReflFromFile'].setvalue(0)
            self.config['GENERAL','fileNames','externalReflectanceFileName'].setvalue('disamar.sim')
            self.config['ADDITIONAL_OUTPUT','additional','refl_Instr_gridSim'].setvalue(1)

    def __call__(self, *extrafiles, **kwargs):
        self.save_config()
        exename = os.path.basename(self.disamar)
        if os.path.dirname(self.output) != "":
            destwd = os.path.abspath(os.path.dirname(self.output))
        else:
            destwd = os.path.abspath(os.getcwd())
        oldcwd = os.getcwd()
        os.chdir(self.tmploc)
        stdoutf = tempfile.TemporaryFile(mode='w+')
        if self._debug:
            print("TEMPORARY_PATH: {0}".format(self.tmploc))
            for item in os.listdir(self.tmploc):
                print("LISTING --- {0}".format(item))
        try:
            starttime = time.time()
            failure = subprocess.call(os.path.join(self.tmploc, exename),
                                      shell=True,
                                      stdout=stdoutf,
                                      stderr=subprocess.STDOUT)
            timelapse = time.time() - starttime

            stdoutf.seek(0)
            for line in stdoutf:
                if 'stopped' in line.lower():
                     raise RuntimeError("Disamar stopped with an error")

            stdoutf.seek(0)
            if failure < 0:
                sys.stderr.write("DISAMAR was terminated by signal {0:d}\n".format(-failure))
                for line in stdoutf:
                    sys.stderr.write(line)
                raise RuntimeError("DISAMAR was terminated by signal {0:d}\n".format(-failure))
            elif failure != 0:
                sys.stderr.write("DISAMAR returned an error ({0:d})\n".format(failure))
                for line in stdoutf:
                    sys.stderr.write(line)
                raise RuntimeError("DISAMAR returned an error ({0:d})\n".format(failure))
        except RuntimeError as e:
            sys.stderr.write("Execution of DISAMAR failed: {0}\nDisamar wrote to the output:\n".format(e))
            stdoutf.seek(0)
            for line in stdoutf:
                sys.stderr.write("\t" + line)
            raise

        os.chdir(destwd)

        attrs = kwargs.copy()
        attrs['uid'] = os.getuid()
        fmt = "{time[0]:04d}{time[1]:02d}{time[2]:02d}T{time[3]:02d}:{time[4]:02d}:{time[5]:02d}"
        attrs['date'] = fmt.format(time=time.gmtime()[0:6])
        attrs['user'] = os.environ['LOGNAME']
        if 'HOSTNAME' in list(os.environ.keys()):
            attrs['host'] = os.environ['HOSTNAME']
        else:
            attrs['host'] = '<Unknown>'

        attrs['time'] = timelapse

        stdoutf.seek(0)
        saw_error = False
        attrs['Error'] = ''
        attrs['ErrorMsg'] = ''

        for line in stdoutf:
            if not self._quiet:
                print(line[0:-1])
            if line.startswith('time for '):
                label, value = line.split('=')
                try:
                    val = float(value)
                except ValueError:
                    val = float('nan')
                attrs["{0}_{1}".format(label[9:-7].replace(' ', '_'), "time")] = val
            if saw_error:
                saw_error = False
                attrs['ErrorMsg'] = line.strip()
            if 'ERROR' in line:
                saw_error = True
                attrs['Error'] = line.strip()
        stdoutf.close()

        if extrafiles is None:
            extrafiles = []
        else:
            extrafiles = list(extrafiles)
        extrafiles.append('disamar.sim')

        self.attributes = attrs

        try:
            if self._debug:
                print("temploc ({0}) listing:".format(self.tmploc))
                for item in os.listdir(self.tmploc):
                    print("\t{0}".format(item))

            translate(os.path.join(self.tmploc, 'disamar.asciiHDF'),
                                self.output,
                                attachment=os.path.join(self.tmploc, 'Config.in'),
                                attributes=attrs, debug=self._debug)
            for name in extrafiles:
                if os.path.exists(os.path.join(self.tmploc, os.path.basename(name))):
                    ext = os.path.splitext(name)[1]
                    shutil.copy(os.path.join(self.tmploc, os.path.basename(name)),
                                os.path.join(destwd,
                                             os.path.splitext(os.path.basename(self.output))[0] + ext))
                else:
                    if os.path.basename(name) != "disamar.sim":
                        print("File {0} not found".format(os.path.basename(name)))
                        print("Temporary location: {0}".format(self.tmploc))
        except (ValueError, IOError) as err:
            print(err)
            for filename in ['disamar.asciiHDF', 'disamar.sim', 'Config.in']:
                ext = os.path.splitext(filename)[1]
                try:
                    shutil.copy(os.path.join(self.tmploc, filename),
                                os.path.join(destwd,
                                  os.path.splitext(os.path.basename(self.output))[0] + ext))
                except:
                    pass
            raise
        os.chdir(oldcwd)

    def make_copy(self):
        '''Make a copy of the DISAMAR executable and support files'''
        disdir = os.path.dirname(os.path.abspath(self.disamar))
        exename = os.path.basename(self.disamar)
        refspecdir = os.path.join(disdir, 'RefSpec')
        scatterdir = os.path.join(disdir, 'expCoefFiles')
        tomsdir = os.path.join(disdir, 'TOMS_v8')

        try:
            os.symlink(refspecdir, os.path.join(self.tmploc, 'RefSpec'))
        except:
            if self._debug:
                sys.write.stderr("Could not create symlink to 'RefSpec'\n")
            raise

        try:
            os.symlink(scatterdir, os.path.join(self.tmploc, 'expCoefFiles'))
        except:
            if self._debug:
                sys.write.stderr("Could not create symlink to 'expCoefFiles'\n"
                                 "Assuming that this version of DISAMAR doesn't need them\n")
        try:
            os.symlink(tomsdir, os.path.join(self.tmploc, 'TOMS_v8'))
        except:
            if self._debug:
                sys.write.stderr("Could not create symlink to ''TOMS_v8''\n"
                                 "Assuming that this version of DISAMAR doesn't need them\n")

        os.symlink(os.path.join(disdir, exename), os.path.join(self.tmploc, exename))

        # shutil.copy(os.path.join(disdir, exename), os.path.join(self.tmploc, exename))




'''
Utilities module
================

  These are some functions to perform some common tasks are
  collected for use in the TropomiDataProcessor collection.

  @author: Maarten Sneep
  @contact: maarten.sneep@knmi.nl
  @organization: KNMI
  @copyright: Maarten Sneep, KNMI
  @license: Creative Commons Attribution-ShareAlike 3.0 Netherlands License.
'''



def classFromString(className, mod=None):
    """Load a module and get the class from a string.

By default the module name is assumed to be the same as the classname.
It is possible to set a separate module name through the mod keyword
parameter.

@param className: The name of the class to load. By default this function
looks for C{Name.Name} (the class C{Name} in the module C{Name}).
This can be changed with the C{mod} parameter.
@type className: C{string}
@param mod: The module name, to override the default.
@type mod: C{string}

@return: The class object that is referenced.
@rtype: A class object.

@note: The function returns C{None} if the module can't be loaded.
ImportError and NameError are handled by the function, other errors are reC{raise}d.
"""
    if mod is None:
        mod = className
    if className == "NoneType":
        cls = None
    else:
        try:
            __import__(mod, globals(), locals(), [], -1)
            cls = sys.modules[mod].__dict__[className]
        except ImportError:
            try:
                cls = eval("{0}".format(className))
            except NameError:
                print('Class "{0}" from modue "{1}"'
                    ' was not found.'.format(className, mod))
                return
            except:
                print('An unanticipated error occurred '
                      'while trying to find Class "{0}"'
                      ' in module "{1}".'.format(className, mod))
                raise
        except:
            print('Module "{0}" was not found, terminating'.format(mod))
            raise
    return cls

def locateResource(name, loc="tbl", isFile=True, mustExist=True, base=None):
    """Locate a resource within the source distribution.

A resoure can be any file or directory in the TropomiDataProcessor
distribution. Before looking in the source distribution, we look in the
current directory. A warning is printed if the file is found there (but
it will be used).

@param name: The name of the file.
@type name: C{string}
@param loc: The subdirectory in which to look for the file.
Defaults to C{tbl} for the C{tbl} sibling of the C{src} directory.
@type loc: C{string}
@param isFile: Is the specified resource a file? Defaults to C{True}.
@type isFile: Bool
@param mustExist: Must the specified resource exist? Defaults to C{True}.
@type mustExist: Bool
@param base: the file which is to be used as the starting point for the search.
    The default is to start from C{__file__}.

@return: The fully qualified path to the specified resource.
@rtype: C{string}

@raise ValueError: If the specified resource does not match the
requirements that are specified with the C{isFile} and C{mustExist}
parameters.
"""
    if mustExist and isFile and os.path.exists(name):
        if os.path.isabs(name):
            path = name
        else:
            path = os.path.realpath(name)
            sys.stderr.write("""Found file "{0}" in the current directory, not searching in "{1}".\n""".format(name, loc))
    else:
        if base is None:
            base = __file__

        path = os.path.join(
                   os.path.dirname(
                       os.path.dirname(
                           os.path.realpath(base))),
                           loc, name)

    if mustExist and (not os.path.exists(path)):
        raise ValueError('File "{0}" not found in the distribution'
            ' (in the "{1}" directory).'.format(name, loc))
    if (mustExist) and (isFile) and (not os.path.isfile(path)):
        raise ValueError('Item "{0}" is not a file.'.format(name))
    if (mustExist) and (not isFile) and (not os.path.isdir(path)):
        raise ValueError('Item "{0}" is a file, expected a directory.'.format(name))

    return path

_basename = 'knmi'

def logname():
    """
    Raise a fake exception, Find the stack frame of the caller so that we
    can traverse to the source filename that started the code.
    From this 'start' code filename the basename is taken without the .py

    If the component found is equal to "pdb", it is assumed the process is
    running inside the python debugger and element 5, instead of 0, of the
    stack is examined to find the component.

    If the component found is equal to "threading", it is assumed the process is
    running inside a thread object and element 1, instead of 0, of the stack is
    examined to find the component.

    Code inspired by (Ian van der Neut).
    """
    global _basename

    parent = os.path.splitext(os.path.basename(wheresdaddy()))[0]
    return '.'.join([_basename, os.path.splitext(os.path.basename(sys.argv[0]))[0], parent])

def whoami():
    return inspect.stack()[1][3]
def whosdaddy():
    return inspect.stack()[2][3]
def whereami():
    return inspect.stack()[1][1]
def wheresdaddy():
    return inspect.stack()[2][1]

def indent_string(s, indent="    "):
    return '\n'.join([indent + c for c in s.split('\n')])

class MyFormatter(logging.Formatter):
    """
    Subclass to change the formatting of the message string.

    Indent every line of the message string.
    """
    def format(self, record):
        """
        Format the specified record as text.

        The record's attribute dictionary is used as the operand to a
        string formatting operation which yields the returned string.
        Before formatting the dictionary, a couple of preparatory steps
        are carried out. The message attribute of the record is computed
        using LogRecord.getMessage(). If the formatting string contains
        "%(asctime)", formatTime() is called to format the event time.
        If there is exception information, it is formatted using
        formatException() and appended to the message.
        """
        record.message = indent_string(record.getMessage())
        if "%(asctime)" in self._fmt:
            record.asctime = self.formatTime(record, self.datefmt)
        s = self._fmt % record.__dict__
        if record.exc_info:
            # Cache the traceback text to avoid converting it multiple times
            # (it's constant anyway)
            if not record.exc_text:
                record.exc_text = self.formatException(record.exc_info)
        if record.exc_text:
            if s[-1:] != "\n":
                s = s + "\n"
            s = "{0}    Exception:\n    {1}".format(s, indent_string(record.exc_text))
        return s

def logger():
    name = logname()

    # create the logging object
    return logging.getLogger(name)

def setup_logging(name=None, handlers=None, levels=None, filenames=None):
    global _basename
    if name is not None and type(name) == type(""):
        _basename = name

    logger = logging.getLogger(_basename)

    logger.setLevel(logging.NOTSET)

    # create formatter
    formatter = MyFormatter(fmt='%(asctime)s %(name)s %(levelname)s %(module)s %(funcName)s %(lineno)d:\n\t%(message)s',
                            datefmt='%Y-%m-%dT%H:%M:%S')
    formatter.converter = time.gmtime

    if handlers is None:
        handlers = [logging.StreamHandler]

    if levels is None:
        levels = ["ERROR"]

    if filenames is None:
        filenames = [None]

    # add all requeste handlers. For now only stream and file handlers are allowed.
    for cls, level, filename in zip(handlers, levels, filenames):
        if cls is logging.StreamHandler:
            h = cls()
        elif cls is logging.FileHandler:
            h = cls(filename, mode='a')
        else:
            continue

        numeric_level = getattr(logging, level.upper(), None)
        h.setLevel(numeric_level)
        h.setFormatter(formatter)
        logger.addHandler(h)

    # capture warnings.
    logging.captureWarnings(True)

def available_memory(kind=None):
    """Returns the RAM (total, used & free) of a linux system"""
    values = [v.split() for v in subprocess.check_output(["/usr/bin/free", "-m"]).split('\n')[1:] if v.split()]

    d = dict([(k[0][0:-1], dict(list(zip(('total', 'used', 'free'), [int(s) for s in k[1:]]))))
              for k
              in [v[0:4]
                  for v
                  in values
                  if v[0] in ('Mem:', 'Swap:')]])
    if kind in ('Mem', 'Swap'):
        return d[kind]
    else:
        return d

'''
A pytables based routine for reading ASCII HDF files,
creating HDF files in the process.

When called as a program it will translate between formats.

@author: Maarten Sneep <maarten.sneep@knmi.nl>
'''

warnings.filterwarnings("ignore", "support for unicode type is very limited, and only works for strings that can be cast as ascii")

# __version = "0.1a"
# __date = "2010-03-11"
# __author = "Maarten Sneep <maarten.sneep@knmi.nl>"

class DatastoreError(Exception):
    '''Error class for all Datastore errors'''
    pass

class DatastoreNode(object):
    '''Base class for the objects in our Datastore hierarchy.

Both groups and leafs inherit from this class, and as a result
all of those can handle attributes.'''
    def __init__(self, name=None, **kwargs):
        '''
        Constructor
        '''
        self._name = name.replace(' ', '_')
        self._parent = None

        if 'attributes' in kwargs.keys():
            self._attributes = kwargs['attributes']
        else:
            self._attributes = {}

        self._debug = 'debug' in kwargs.keys() and kwargs['debug']

    def name(self):
        '''Get the name of the node'''
        return self._name

    def set_name(self, name):
        '''Set the name of the node.'''
        self._name = name

    def set_parent(self, parent):
        '''Set a reference to the parent'''
        self._parent = parent

    def parent(self):
        '''Return a reference to the parent node'''
        return self._parent

    def root(self):
        '''Return a reference to the root object'''
        return self.parent().root()

    def attributes(self):
        '''Return all attributes'''
        return self._attributes

    def set_attributes(self, attrs):
        '''Set attributes to attrs. Replaces all existing attributes'''
        self._attributes = attrs

    def set_attribute(self, key, val):
        '''set a single attribute value'''
        self._attributes[key] = val

    def append_attributes(self, **kwargs):
        '''Add keyword arguments to the attributes'''
        self._attributes.update(kwargs)

    def replace_attributes(self, **kwargs):
        '''Replace existing attributes with keyword values'''
        self.set_attributes(kwargs)

    def attribute(self, key):
        '''Return the value of an attribute'''
        return self._attributes[key]

    def leaf(self):
        '''Determine if the instance is a leaf object'''
        return isinstance(self, DatastoreLeaf)


class DatastoreGroup(DatastoreNode):
    '''Base class for Datastore objects that can contain other groups or leafs'''
    def __init__(self, name=None, **kwargs):
        super(DatastoreGroup, self).__init__(name=name, **kwargs)
        self._children = []
        if 'children' in kwargs.keys():
            for child in kwargs['children']:
                self[child.name()] = child

    def __iter__(self):
        return self.next()

    def next(self):
        '''Iterate over all children '''
        for child in self._children:
            yield self.__dict__[child]

    def children(self):
        '''Return a dict with children of this group'''
        children = {}
        for name in self._children:
            children[name] = self.__dict__[name]
        return children

    def add_child(self, child):
        '''Add a child (group or leaf) to the group'''
        self[child.name()] = child

    def items(self):
        '''Return all child items'''
        return self._children[:]

    def keys(self):
        '''return the available child names'''
        return self._children

    def __len__(self):
        return len(self._children)

    def __getitem__(self, key):
        print(key)
        return self.__dict__[key]

    def __setitem__(self, key, val):
        if not issubclass(val.__class__, DatastoreNode):
            raise DatastoreError('The value must be a subclass of DatastoreNode, ' +
                             'instead it is of class {0}'.format(val.__class__.__name__))
        if key != val.name():
            raise KeyError('The key must be equal to the name of the child')
        self.__dict__[key] = val
        self.__dict__[key].set_parent(self)
        self._children.append(key)

    def __contains__(self, item):
        return 1 if item in self._children else 0

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        else:
            raise AttributeError('Attribute {0} not found'.format(name))

    def __setattr__(self, name, val):
        if name in self.__dict__ or name[0] == '_':
            self.__dict__[name] = val
        else:
            self[name] = val

    def __delattr__(self, name):
        if name in self._children:
            del self.__dict__[name]
            self._children.remove(name)
        elif name in self.__dict__:
            del self.__dict__[name]
        else:
            raise AttributeError('Attribute {0} not found'.format(name))


class Datastore(DatastoreGroup):  # IGNORE:R0904
    '''The "root" object. Here the filehandle is maintained'''
    def __init__(self, name=None, **kwargs):
        '''Initialize a datagroup.

The 'name' keyword is used to give a name to the whole dataset. Defaults to 'root'.

If the read keyword is set to a filename, then the data in the group will
be initialized from this file.

If the write keyword is set to a filename, then the data (provided with the 'children'
keyword) are written to that file.

The group is initialized with its children. If the 'children' keyword is set to a
list of DatastoreGroup or DatastoreLeaf objects, then the data object is initialized
with it.

Attributes for the group can be supplied as a dictionary to the 'attributes' keyword.
        '''
        if name is None:
            name = 'root'
        super(Datastore, self).__init__(name=name, **kwargs)
        self._ishdf = None
        self._filename = None
        self._filehandle = None

        keys = kwargs.keys()
        if 'read' in keys:
            self.read(kwargs['read'])
        elif 'write' in keys:
            self.write(kwargs['write'])

    def read(self, fname):
        '''Read data from file'''
        self._filename = fname
        if os.path.exists(self._filename):
            self._ishdf = tables.is_hdf5_file(self._filename)
        else:
            raise DatastoreError("File {0} not found".format(self._filename))

        if self._ishdf:
            if self._debug:
                print('Reading from HDF-5', file=sys.stderr)
            self._filehandle = tables.open_file(self._filename, mode='r')
            self.read_hdf()
        else:
            if self._debug:
                print('Reading from ASCII', file=sys.stderr)
            self._filehandle = open(self._filename, 'r')
            self.read_ascii()
        self._filehandle.close()
        self._filehandle = None

    def write(self, fname):
        '''write data to file'''
        self._filename = fname
        if os.path.exists(self._filename):
            self._ishdf = tables.is_hdf5_file(self._filename)
        else:
            self._ishdf = (os.path.splitext(self._filename)[1] in
                           ['.h5', '.hdf5', '.hdf', '.he5'])

        if self._ishdf:
            if self._debug:
                print('Writing to HDF-5', file=sys.stderr)
            filters = tables.Filters(complevel=6, complib='zlib', fletcher32=True)
            self._filehandle = tables.open_file(self._filename, mode='w',
                                               title=self.name(), filters=filters)
            self.write_hdf()
        else:
            if self._debug:
                print('Writing to ASCII', file=sys.stderr)
            self._filehandle = open(self._filename, 'w')
            self.write_ascii()
        self._filehandle.close()
        self._filehandle = None

    def attach_file(self, file_to_attach, attributes=None): #IGNORE:R0912
        '''Attach a file to the hdf-5 file.

The file is attached to the '/attachments' group.'''
        if not self._ishdf:
            raise DatastoreError('Files can only be attached to HDF-5 files')
        if not os.path.exists(file_to_attach):
            raise ValueError('File not found: {0}'.format(file_to_attach))
        if not os.path.exists(self._filename):
            raise DatastoreError('File must exist before ' +
                                 'attaching: {0}'.format(self._filename))

        if self._debug:
            print('Attaching file "{0}" to HDF-5'.format(file_to_attach), file=sys.stderr)

        maxlen = 65535
        filters = tables.Filters(complevel=9, complib='zlib', fletcher32=True)
        atom = tables.StringAtom(itemsize=maxlen, shape=(), dflt='')

        self._filehandle = tables.open_file(self._filename, mode='a')
        try:
            gnode = self._filehandle.get_node('/attachments')
        except  tables.NoSuchNodeError:
            gnode = self._filehandle.create_group('/', 'attachments', 'file attachments')

        with open(file_to_attach, 'r') as fp:
            text = fp.read()

        numattach = 1 + (len(text) // maxlen)
        stop = len(text)

        for ii in range(numattach):
            name = os.path.basename(file_to_attach)
            name = name.replace('.', '_').replace(' ', '_')
            name = name.replace('-', '_').replace('+', '_')
            if numattach > 1:
                name = '{0}_{1:02d}'.format(name, ii + 1)
            carray = self._filehandle.create_carray(gnode, name, atom, (1,), filters=filters)

            # split on whole lines only.
            start = 0 if ii == 0 else stop
            if ii < numattach - 1:
                stop = (ii + 1) * maxlen
                while text[stop - 1] != '\n' and stop > start:
                    stop -= 1
            else:
                stop = len(text)

            carray[:] = bytes(text[start:stop], 'ASCII')
            carray.flush()
            if numattach > 1:
                carray._v_attrs['Datastore_part'] = ii + 1

            if attributes:
                for k, v in attributes.items():
                    if type(v) == type(u''):
                        v = str(v)
                    elif type(v) == type(True):
                        v = 'true' if v else 'false'
                    elif type(v) == type([]):
                        if type(v[0]) == type(u''):
                            v = [str(e) for e in v]
                        v = np.asarray(v)
                    carray._v_attrs[k] = v



        self._filehandle.close()
        self._filehandle = None

    def get_attachment(self, attachment_name):
        '''Get an attachment from the hdf-5 file'''
        if not self._ishdf:
            raise DatastoreError('Files can only be attached to HDF files')
        if not os.path.exists(self._filename):
            raise DatastoreError('File must exist before reading ' +
                                 'attachment: {0}'.format(self._filename))

        self._filehandle = tables.open_file(self._filename, mode='r')

        try:
            gnode = self._filehandle.get_node('/attachments')
        except  tables.NoSuchNodeError:
            raise DatastoreError('No attachments were found')

        name = os.path.basename(attachment_name)
        name = name.replace('.', '_').replace(' ', '_')
        name = name.replace('-', '_').replace('+', '_')

        attributes = {}
        text = ''

        for item in gnode:
            if item._v_name.startswith(name):
                if not attributes:
                    for k in item.attrs._f_list():
                        if k != 'Datastore_part':
                            attributes[k] = item.attrs[k]
                text = text + item.read().tostring().strip('\x00')

        self._filehandle.close()
        self._filehandle = None

        return (text, attributes)

    def read_ascii(self):
        '''Read a data-object from ascii file'''
        ascii_read_object = AHDFr(fp=self._filehandle, debug=self._debug)
        dstore = ascii_read_object()
        _translate_old_datastore(dstore, self)

    def write_ascii(self):
        '''Write a data-object to ascii file'''
        _write_ascii(self._filehandle, self)

    def read_hdf(self):
        '''Read a data-object from hdf-5/PyTables file'''
        _translate_pytables(self._filehandle.root, self)

    def write_hdf(self):
        '''Write a data-object to hdf-5/PyTables file'''
        _write_hdf5(self._filehandle, self)

    def __del__(self):
        if self._filehandle is not None:
            self._filehandle.close()

    def root(self):
        '''Return a reference to the root object (endpoint of call chain)'''
        return self


class DatastoreLeaf(DatastoreNode):
    '''Leaf object in our datastore hierarchy'''
    def __init__(self, name=None, **kwargs):
        '''Leaf object.

The name of the object is given in the 'name' keyword.

Data

Attributes for the leaf can be supplied as a dictionary to the 'attributes' keyword.
       '''
        super(DatastoreLeaf, self).__init__(name=name, **kwargs)

        self._data = None
        if 'data' in kwargs.keys():
            self._data = kwargs['data']

    def __getitem__(self, key):
        return self._data.__getitem__(key)

    def __setitem__(self, key, val):
        self._data.__setitem__(key, val)

    def __setattr__(self, name, val):
        if name in self.__dict__ or name[0] == '_':
            self.__dict__[name] = val
        else:
            self._data.__setattr__(name, val)

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        else:
            return self._data.__getattribute__(name)

    def __delattr__(self, name):
        self._data.__delattr__(name)

    def data(self):
        '''Return the data object (numpy array)'''
        return self._data

    def set_data(self, newdata):
        '''Change the data array'''
        self._data = newdata

    def keys(self): #IGNORE:R0201
        '''Instead of letting __getattrs__ fail,
we raise a warning here to indicate that we are a leaf'''
        raise UserWarning("Arrived at leaf of object tree")

def read(filename, debug=False):
    '''Read a structured ASCII or HDF-5 input file'''
    root = Datastore(name='root', read=filename, debug=debug)
    return root

def write(filename, obj):
    '''Write a Datastore object to file'''
    obj.write(filename)

def translate(infile, outfile, attachment=None, attributes=None, debug=False):
    '''Translate a file (asciiHDF or HDF-5) into another format (asciiHDF or HDF-5).

Types are determined from file or filename.'''
    d = read(infile, debug=debug)
    write(outfile, d)
    if attachment is not None and d._ishdf:  # IGNORE:W0212
        if type(attachment) == type(""):
            attachment = [attachment]
        for item in attachment:
            if os.path.exists(item) and os.path.isfile(item):
                d.attach_file(item, attributes=attributes)
            else:
                print("Item '{0}' could not be attached because it doesn't exist or is a directory.".format(item))

def _translate_pytables(node, dest):
    '''Translate a pytables node into a hierarchy of DatastoreNodes'''
    for name in node._v_attrs._f_list(): # IGNORE:W0212
        dest.set_attribute(name, node._v_attrs[name]) # IGNORE:W0212
    for item in node:
        if issubclass(item.__class__, tables.Leaf):
            attrs = {}
            for name in item._v_attrs._f_list(): # IGNORE:W0212
                attrs[name] = item._v_attrs[name] # IGNORE:W0212
            child = DatastoreLeaf(name=item._v_name, # IGNORE:W0212
                                  data=item.read(),
                                  attributes=attrs)
            dest[child.name()] = child
        else:
            child = DatastoreGroup(name=item._v_name) # IGNORE:W0212
            dest[child.name()] = _translate_pytables(item, child)
    return dest

def _translate_old_datastore(dstore, dest):
    '''Translate an old Datastore object into a hierarchy of DatastoreNodes'''
    for key, val in dstore.items():
        dest.set_attribute(key, val)
    for item in dstore.data():
        if item.leaf():
            child = DatastoreLeaf(name=item.name(),
                                  data=item.data(),
                                  attributes=item.copy())
            dest[child.name()] = child
        else:
            child = DatastoreGroup(name=item.name())
            dest[child.name()] = _translate_old_datastore(item, child)
    return dest

class _write_hdf5(object):
    '''Write data to an HDF-5 file'''
    def __init__(self, fp=None, datastore=None):
        self.datastore = datastore
        self.fp = fp
        self.write_group(self.datastore, fp.root)

    def write_group(self, datastore, parent):
        '''Write datastore as a group.'''
        if not datastore.leaf():
            if datastore.name() == 'root':
                group = self.fp.root
            else:
                group = self.fp.create_group(parent, datastore.name(), datastore.name())
            self.write_attributes(datastore, group)
            self.write_data(datastore, group)

    def write_attributes(self, datastore, parent):  #IGNORE:R0201
        '''Write the key-value pairs in datastore.attributes() to an attributes block'''
        attrs = datastore.attributes()
        for k, v in attrs.items():
            # make sure that attribute names do not contain spaces.
            if ' ' in k:
                k = k.replace(' ', '_')

            # do some type cleanup
            if type(v) == type(u''):
                v = str(v)
            elif type(v) == type(True):
                v = 'true' if v else 'false'
            elif type(v) == type([]):
                if type(v[0]) == type(u''):
                    v = [str(e) for e in v]
                v = np.asarray(v)

            parent._v_attrs[k] = v  #IGNORE:W0212

    def write_data(self, datastore, parent):
        '''Write the data in datastore.

The datastore can be a group (in which case we dispatch to write_group)
or a leaf, and then we write the array.
'''
        if datastore.leaf():
            self.write_array(datastore, parent)
        else:
            for item in datastore:
                if item.leaf():
                    self.write_array(item, parent)
                else:
                    self.write_group(item, parent)

    def write_array(self, datastore, parent):
        '''Write the array to the file'''
        atom = tables.Atom.from_dtype(datastore.data().dtype)
        shape = datastore.data().shape

        if 0 in shape:
            return

        try:
            ca = self.fp.create_carray(parent, datastore.name(), atom, shape)
        except tables.exceptions.NodeError as err:
            print("Warning: {0}".format(err))
            ca = self.fp.get_node(parent, datastore.name())

        ca[...] = datastore.data()

        self.write_attributes(datastore, ca)

class _write_ascii(object):
    '''Write data to an ASCII file'''
    def __init__(self, fp=None, datastore=None):
        self.datastore = datastore
        self.fp = fp
        self.write_group(self.datastore)

    def write_group(self, datastore):
        '''Write datastore as a group.'''
        if not datastore.leaf():
            self.fp.write("BeginGroup({0})\n".format(datastore.name()))
            self.write_attributes(datastore)
            self.write_data(datastore)
            self.fp.write("EndGroup\n")

    def write_attributes(self, datastore):
        '''Write the key-value pairs in datastore to an attributes block'''
        attrs = datastore.attributes()
        if len(attrs) > 0:
            self.fp.write("BeginAttributes\n")
            for key, val in attrs.items():
                self.fp.write("{key} = {value}\n".format(key=key,
                              value=json.dumps(val)))
            self.fp.write("EndAttributes\n")

    def write_data(self, datastore):
        '''Write the data in datastore.

The datastore can be a group (in which case we dispatch to write_group)
or a leaf, and then we write the array.
'''
        if datastore.leaf():
            self.write_array(datastore)
        else:
            for item in datastore:
                if item.leaf():
                    self.write_array(item)
                else:
                    self.write_group(item)

    def write_array(self, data):
        '''Write the array to the file.

The attributes in the datastore are added to the array block. The order (only C),
number of dimensions, and dimensions sizes are determined from the array object.
The numpy tofile method is used to add the data itself.
'''
        if (str(data.data().dtype)[0:2].lower() != '|s'):
            typedescr = str(data.data().dtype)
        else:
            typedescr = 'String'

        self.fp.write("BeginArray({name}, {type})\n".format(name=data.name(),
                                                            type=typedescr))
        self.write_attributes(data)
        self.fp.write("Order = C\n")
        self.fp.write("NumDimensions = {ndims}\n".format(ndims=len(data.data().shape)))
        self.fp.write("Size = {0}\n".format(", ".join([str(d) for d in data.data().shape])))
        if typedescr != 'String':
            data.data().tofile(self.fp, sep=' ')
        else:
            data.data().tofile(self.fp, sep='\n')
        self.fp.write("\nEndArray\n")


'''
This module can be used for data storage, and reading and
writing to a simple ascii format.
'''



# __date = "2010-02-26"
# __version = "1.0"
# __author = "Maarten Sneep"


class _Datastore(dict): # IGNORE:R0904
    '''Data storage object to hold data (groups, arrays) and attributes.

The data is stored as some object:
  - list of _Datastore objects for a group.
  - A numpy array for a leaf object.
the attributes are stored in the inherited dict.
    '''

    def __init__(self, name=None, data=None, attrs=None):
        """Initializer

Copy the data object, and put the attributes in the inherited
dictionary. Since we inherit from dict, all normal methods to handle
key-value pairs are available.
"""
        if not attrs:
            attrs = {}
        dict.__init__(self, attrs)
        self._name = name
        self._data = data
        self._leaf = None
        self.set_leaf()

    def set_leaf(self):
        '''Boolean to indicate whether this is a leaf or a group.

This is automatically determined from the type of data.
- List => group
- Array => leaf
'''
        if self.data is None:
            self._leaf = None
        else:
            if type(self._data) == type([]):
                self._leaf = False
            else:
                self._leaf = True

    def leaf(self):
        '''Is this a leaf or a group?'''
        self.set_leaf()
        return self._leaf

    def attributes(self):
        """Return attributes for the data object (i.e. return self)"""
        return self

    def set_attributes(self, attrs):
        """Replace all attributes"""
        self.clear()
        self.update(attrs)

    def data(self):
        """Return the data object (numpy array or list)"""
        return self._data

    def set_data(self, data):
        """Replace the data object (numpy array or list)"""
        if self._data:
            del self._data
        self._data = data
        self.set_leaf()

    def item(self, aname, avalue):
        """Return the _Datastore objects where attribute 'aname' is set to 'avalue'.

Return value is a tuple of found objects.
"""
        rval = []
        if aname in list(self.keys()) and self[aname] == avalue:
            rval.append(self)

        for obj in self.data():
            if obj.leaf() and aname in list(obj.keys()) and obj[aname] == avalue:
                rval.append(obj)
            else:
                rval.extend(obj.item(aname, avalue))

        return tuple(rval)

    def append(self, data):
        """Add a _Datastore object to the group"""
        if self.leaf():
            raise ValueError("Cannot append to a leaf object")
        if type(data) == type([]):
            self._data.extend(data)
        else:
            self._data.append(data)

    def name(self):
        return self._name

class AHDFr(object):
    """Object to read ASCII HDF representation"""
    def __init__(self, fname=None, fp=None, root=None, debug=False):
        self.re = {}
        self.re['begingroup'] = re.compile(r"^\s*BeginGroup\(([^)]+)\)")
        self.re['endgroup'] = re.compile(r"^\s*EndGroup")
        self.re['beginattrs'] = re.compile(r"^\s*BeginAttributes")
        self.re['attr'] = re.compile(r"^\s*([a-zA-Z][^=]*?)=(.+?)$")
        self.re['endattrs'] = re.compile(r"^\s*EndAttributes")
        self.re['beginarray'] = re.compile(r"^\s*BeginArray\(([^,]+),\s*([^)]+)\)")
        self.re['endarray'] = re.compile(r"^\s*EndArray")
        self.re['comment'] = re.compile(r"^\s*#")
        self.re['order'] = re.compile(r"^\s*Order\s*=\s*(\w+)")
        self.re['ndim'] = re.compile(r"^\s*NumDimensions\s*=\s*(\d+)")
        self.size = lambda n: re.compile(r"^\s*Size\s*=\s*{0}".format(
                                         r"\s*[, ]\s*".join([r"(\d+)"] * n)))
        self.fp = None

        if root is None:
            self.root = _Datastore(name="/", data=[], attrs={})
        else:
            self.root = root

        if fname is not None:
            self.fp = open(fname, "r")
        if fp is not None and self.fp is None:
            self.fp = fp
        self._debug = debug
        self.props = {'continue_scan':True}

    def __del__(self):
        if self.root.name() == "/":
            self.fp.close()

    def __call__(self):
        '''Read a group from the file.

The search routine will reset the continue_scan property when the group is closed.
'''
        self.props['continue_scan'] = True
        while self.props['continue_scan']:
            line = self.fp.readline()
            if not line:
                self.props['continue_scan'] = False
            else:
                for key in ['begingroup', 'endgroup',
                            'beginattrs', 'endattrs',
                            'beginarray', 'endarray',
                            'comment']:
                    search = self.re[key].search(line)
                    if search:
                        if self._debug:
                            print('Found key: {0}'.format(key), file=sys.stderr)
                        f = eval("self." + key)
                        f(search)
                        break

        return self.root

    def begingroup(self, search):
        '''Handle the BeginGroup statement in the data file.

If the name is not "/", then a new AHDF object is created,
and the group is read in this new object.
'''
        self.props['group'] = True

        name = search.group(1)
        if self._debug:
            print('Found group: {0}'.format(name), file=sys.stderr)
        if not (name == "/" or name == ""):
            basename = os.path.basename(name)
            obj = AHDFr(fp=self.fp, root=_Datastore(name=basename, data=list()), debug=self._debug)
            self.root.append(obj())

    def endgroup(self, dummy):
        '''Handle the EndGroup directive

Set continue_scan property to False to end scanning the current group.
'''
        self.props['continue_scan'] = False

    def beginattrs(self, dummy, array=False):
        '''Handle the BeginAttributes directive'''
        attrs = {}
        if self._debug:
            print('Reading {0} attributes'.format('array' if array else 'group'), file=sys.stderr)
        while True:
            before = self.fp.tell()
            line = self.fp.readline()
            if len(line) == 0 or line == '\n' or self.re["endattrs"].search(line):
                self.fp.seek(before)
                break
            else:
                search = self.re["attr"].search(line)
                if search is None and self._debug:
                    print('Attribute search failed: {0}'.format(line[0:-1]))
                try:
                    attrs[search.group(1).strip()] = json.loads(search.group(2).strip())
                except:
                    print('Encountered error while parsing')
                    print(line)
                    raise

        if array:
            self.root.data()[-1].set_attributes(attrs)
        else:
            self.root.update(attrs)

    def attr(self, dummy):# IGNORE:R0201
        '''read attribute. Completely handled in self.beginattrs()'''
        raise ValueError("Unexpected call")

    def endattrs(self, dummy):
        '''End the attributes directive (No-op)'''
        pass

    def beginarray(self, search):
        '''Start reading a new array.

This method handles the array header, and records the name and type of the array.
The data itself is read with the readarraydata() method.
'''
        self.props['name'] = search.group(1).strip()
        self.props['type'] = search.group(2).strip()
        self.props['group'] = False
        if self._debug:
            print('Start reading array "{0}"'.format(search.group(1).strip()), file=sys.stderr)

        self.root.append(_Datastore(name=self.props['name'],
                                    data=np.zeros((0,), dtype='float64')))
        self.readarraydata()

    def endarray(self, dummy):
        '''Remove the properties to prepare the object for the next array'''
        pass

    def readarraydata(self):
        '''Read the data.

Start with the header (Order, NumDimensions, and Size keywords),
and finally read the data itself.
'''

        if self.props['group']:
            if self._debug:
                print('reading array aborted', file=sys.stderr)
            return

        # read array attributes
        line = self.fp.readline()
        search = self.re['beginattrs'].search(line)
        if search:
            self.beginattrs(line, array=True)
            line = self.fp.readline()
            self.endattrs(line)
            line = self.fp.readline()

        # read the order (Fortran or C)
        search = self.re["order"].search(line)
        if not search:
            raise ValueError("Unexpected key:\n\t'{0}'\n\tExpected 'Order'".format(line))
        self.props['reverse'] = (search.group(1).lower() == "fortran")
        if self._debug:
            print('Array order "{0}"'.format(search.group(1)), file=sys.stderr)

        # read the number of dimensions
        line = self.fp.readline()
        search = self.re["ndim"].search(line)
        if not search:
            raise ValueError("Unexpected key:\n\t'{0}'"
                             "\n\tExpected 'NumDimensions'".format(line))
        self.props['ndim'] = int(search.group(1))
        if self._debug:
            print('Number of dimensions "{0}"'.format(int(search.group(1))), file=sys.stderr)

        # read the dimension sizes
        line = self.fp.readline()
        search = self.size(self.props['ndim']).search(line)
        if not search:
            raise ValueError("Unexpected key:\n\t'{0}'\n\tExpected 'Size'".format(line))
        self.props['dimensions'] = np.asarray([int(v) for v in search.groups()])

        dims = self.props['dimensions']
        if self._debug:
            print('Dimensions "{0}"'.format(dims.tolist()), file=sys.stderr)

        # String arrays cannot be read directly by numpy,
        # we have to do this with an intermediate list.
        if self.props['type'].lower() == "string":
            array = []
            for dummy in range(dims.prod()):
                array.append(self.fp.readline().strip())

            if self.props['ndim'] > 1:
                if self.props['reverse']:
                    array = np.asarray(array).reshape(dims[::-1])
                    array.transpose()
                else:
                    array = np.asarray(array).reshape(dims)
            else:
                array = np.asarray(array)
        else:
            array = np.fromfile(file=self.fp, dtype=np.dtype(self.props['type']),
                                    count=dims.prod(),
                                    sep=' ')

            if len(array) != dims.prod():
                if self._debug:
                    print('Array {0} has incorrect size.'.format(self.props['name']), file=sys.stderr)
                    oldsize = len(array)
                    array = np.resize(array, (dims.prod(),))
                    array[oldsize:] = 0.0
                    # finish reading the line
                    line = self.fp.readline()
                    print('Rest of the line:\n"{0}"'.format(line), file=sys.stderr)
                    items = line.split()
                    cnt = -1
                    for item in items[::-1]:
                        try:
                            array[cnt] = float(item)
                        except ValueError:
                            array[cnt] = np.nan
                        cnt -= 1
                else:
                    raise ValueError('Unexpected array size for {0}'.format(self.props['name']))

            if self.props['ndim'] > 1:
                if self.props['reverse']:
                    array = np.reshape(array, dims[::-1])
                    array.transpose()
                else:
                    array = np.reshape(array, dims)

        self.root.data()[-1].set_data(array)

    def comment(self, dummy):
        '''Simply ignore comments'''
        pass

    def order(self, dummy): # IGNORE:R0201
        '''This method should not be called directly'''
        raise ValueError("Unexpected call")

    def ndim(self, dummy): # IGNORE:R0201
        '''This method should not be called directly'''
        raise ValueError("Unexpected call")

