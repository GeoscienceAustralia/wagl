#!/usr/bin/env python

"""XML Metadata module

Author: Alex Ip (alex.ip@ga.gov.au)
"""

import xml.dom.minidom
from StringIO import StringIO
import logging, os, re

from ULA3.utils import unicode_to_ascii
from . import Metadata

logger = logging.getLogger('root.' + __name__)

class XMLMetadata(Metadata):
    """Subclass of Metadata to manage XML data
    """
    # Class variable holding metadata type string
    _metadata_type_id = 'XML'
    _filename_pattern = '.*\.xml' # Default RegEx for finding metadata file.

    def __init__(self, source = None, uses_attributes=False):
        """Instantiates XMLMetadata object. Overrides Metadata method
        """
        self._uses_attributes = uses_attributes # Boolean flag indicating whether values are stored as tag attributes
        self.processing_instruction = {} # Dict containing processing instruction name and value
        self.document_attributes = {} # Dict containing any attributes when not self._uses_attributes
        Metadata.__init__(self, source); # Call inherited constructor

    #===========================================================================
    # # Attempted fix for whitespace issue using PyXML. Caused problems with SceneDataset config XML
    # def toprettyxml_fixed (self, node, encoding='utf-8'):
    #     """
    #     Creates well-formatted XML without whitespaces
    #     """
    #     tmpStream = StringIO()
    #     PrettyPrint(node, stream=tmpStream, encoding=encoding)
    #     return tmpStream.getvalue()
    #===========================================================================

    def _populate_dict_from_node(self, node, tree_dict, level=0):
        """Private recursive function to populate a nested dict from DOM tree or element node
        Exposed to allow unit testing using a DOM tree constructed from a string
        Arguments:
            node: xml.dom.Node object to traverse
            tree_dict: nested dict structure to hold result
        """
        # Traverse all non-text nodes
        for child_node in [x for x in node.childNodes if x.nodeType == xml.dom.minidom.Node.ELEMENT_NODE]:
            nodeName = child_node.nodeName
            if nodeName:
                nodeName = unicode_to_ascii(nodeName)

                logger.debug('%sDOM Node name = %s, Node type = %s, Child nodes = %s, Attributes = %s',
                             '  ' * level, nodeName, child_node.nodeType, child_node.childNodes, child_node.attributes)

                subtree_dict = {}
                if child_node.childNodes: # Recursive call to check for non-text child nodes
                    self._populate_dict_from_node(child_node, subtree_dict, level + 1)

                logger.debug('%s  subtree_dict = %s', '  ' * level, subtree_dict)
                if child_node.attributes:
                    logger.debug('%s  Child node attribute count = %s', '  ' * level, len(child_node.attributes))

                if subtree_dict: # Not a leaf node - sub-nodes found
                    tree_dict[nodeName] = subtree_dict

                elif child_node.attributes: # Leaf node - values held in attributes
                    self._uses_attributes = True # Remember that attributes are being used for this file

                    subtree_dict = {}
                    tree_dict[nodeName] = subtree_dict
                    level += 1
                    for attr_index in range(len(child_node.attributes)):
                        attribute = child_node.attributes.item(attr_index)
                        logger.debug('%s  Attribute: %s = %s', '  ' * level, attribute.name, attribute.value)
                        subtree_dict[unicode_to_ascii(attribute.name)] = unicode_to_ascii(attribute.value)
                    level -= 1

                elif child_node.childNodes and child_node.childNodes[0].nodeType == xml.dom.minidom.Node.TEXT_NODE: # Leaf node - value held in child text node
                    # Take value of first text child node
                    tree_dict[nodeName] = unicode_to_ascii(child_node.childNodes[0].nodeValue)
                    logger.debug('%s  Node value = %s from text child node', '  ' * level, tree_dict[nodeName])
                elif not child_node.childNodes: # Empty leaf node
                    tree_dict[nodeName] = ''

    def _populate_node_from_dict(self, tree_dict, node, uses_attributes, owner_document=None, level=0):
        """Private recursive function to populate a nested dict from DOM tree or element node
        Exposed to allow unit testing using a DOM tree constructed from a string
        Arguments:
            tree_dict: nested dict structure to traverse
            node: xml.dom.Node object to hold result
            uses_attributes: Boolean flag indicating whether to write values to tag attributes
        """
        owner_document = owner_document or node

        for node_name in sorted(tree_dict.keys()):
            child_item = tree_dict[node_name]
            assert child_item is not None, node_name + ' node is empty - must hold either a string or subtree dict'
            if type(child_item) == dict: # Subtree - descend to next level
                logger.debug('%sElement Node %s', '  ' * level, node_name)
                child_node = xml.dom.minidom.Element(node_name)
                child_node.ownerDocument = owner_document
                node.appendChild(child_node)

                self._populate_node_from_dict(child_item, child_node, uses_attributes, owner_document, level + 1)

            else: # Leaf node - store node value
                if child_item is None:
                    child_item = ''
                assert type(child_item) == str, node_name + ' node is not a string'
                if uses_attributes: # Store value in attribute
                    logger.debug('%sAttribute for %s = %s', '  ' * level, node_name, repr(child_item))
                    node.setAttribute(node_name, child_item)
                else: # Store value in child text node
                    logger.debug('%sText Child Node for %s = %s', '  ' * level, node_name, repr(child_item))

                    child_node = xml.dom.minidom.Element(node_name)
                    child_node.ownerDocument = owner_document
                    node.appendChild(child_node)

                    # Only add text node if value is non-empty
                    if child_item:
                        text_node = xml.dom.minidom.Text()
                        text_node.ownerDocument = owner_document
                        text_node.nodeValue = child_item
                        child_node.appendChild(text_node)


    def read_file(self, filename=None):
        """Function to parse an XML metadata file and store the results in self._metadata_dict
        Argument:
            filename: XML Metadata file to be parsed and stored
        Returns:
            nested dict containing metadata
        """
        logger.debug('read_file(%s) called', filename)

        filename = filename or self._filename
        assert filename, 'Filename must be specified'

        logger.debug('Parsing XML file %s', filename)

        # Open XML document using minidom parser
        dom_tree = xml.dom.minidom.parse(filename)

        # Remember any processing instruction node
        for node in dom_tree.childNodes:
            if node.nodeType == xml.dom.minidom.Node.PROCESSING_INSTRUCTION_NODE:
                processing_instruction_node = node
                logger.debug('Processing Instruction Node found: Name = %s, Value = %s', processing_instruction_node.nodeName, processing_instruction_node.nodeValue)
                self.processing_instruction['name'] = processing_instruction_node.nodeName
                self.processing_instruction['value'] = processing_instruction_node.nodeValue
            elif node.nodeType == xml.dom.minidom.Node.ELEMENT_NODE: # Root node has attributes
                for attr_index in range(len(node.attributes)):
                    attribute = node.attributes.item(attr_index)
                    logger.debug('Document Attribute: %s = %s', attribute.name, attribute.value)
                    self.document_attributes[unicode_to_ascii(attribute.name)] = unicode_to_ascii(attribute.value)

        # Create nested dict from DOM tree
        self._populate_dict_from_node(dom_tree, self._metadata_dict)
        self._filename = filename

        return self._metadata_dict

    def write_file(self, filename=None, uses_attributes=None, save_backup=False):
        """Function write the metadata contained in self._metadata_dict to an XML file
        Argument:
            filename: Metadata file to be written
            uses_attributes: Boolean flag indicating whether to write values to tag attributes
        """
        logger.debug('write_file(%s) called', filename)

        filename = filename or self._filename
        assert filename, 'Filename must be specified'

        # Allow values to be stored as attributes
        if uses_attributes is None:
            uses_attributes = self._uses_attributes

        if save_backup and os.path.exists(filename + '.bck'):
            os.remove(filename + '.bck')

        if os.path.exists(filename):
            if save_backup:
                os.rename(filename, filename + '.bck')
            else:
                os.remove(filename)

        # Open XML document
        try:
            outfile = open(filename, 'w')
            assert outfile is not None, 'Unable to open XML file ' + filename + ' for writing'

            logger.debug('Writing XML file %s', filename)

            dom_tree = xml.dom.minidom.Document()

            # Write any processing instruction node
            if self.processing_instruction:
                processing_instruction_node = xml.dom.minidom.ProcessingInstruction(self.processing_instruction['name'], self.processing_instruction['value'])
                processing_instruction_node.ownerDocument = dom_tree
                dom_tree.appendChild(processing_instruction_node)

            # Open XML document using minidom parser
            self._populate_node_from_dict(self._metadata_dict, dom_tree, uses_attributes)

            # Set root node attributes if required
            if self.document_attributes:
                root_node = [node for node in dom_tree.childNodes if node.nodeType == xml.dom.minidom.Node.ELEMENT_NODE][0]
                for attribute_name in self.document_attributes.keys():
                    root_node.setAttribute(attribute_name, self.document_attributes[attribute_name])

#            outfile.write(dom_tree.toxml(encoding='utf-8'))
#            outfile.write(dom_tree.toprettyxml(encoding='utf-8'))
#            outfile.write(self.toprettyxml_fixed(node, encoding='utf-8')) # PyXML required

            #===================================================================
            # # Strip all tabs and EOLs from around values
            # outfile.write(re.sub('(\<\w*[^/]\>)\n(\t+\n)*(\t*)([^<>\n]*)\n\t*\n*(\t+)(\</\w+\>)',
            #                      '\\1\\4\\6',
            #                      dom_tree.toprettyxml(encoding='utf-8')
            #                      )
            #               )
            #===================================================================

            # Strip all tabs and EOLs from around values, remove all empty lines
            outfile.write(re.sub('\>(\s+)(\n\t*)\<',
                                 '>\\2<',
                                 re.sub('(\<\w*[^/]\>)\n(\t*\n)*(\t*)([^<>\n]*)\n\t*\n*(\t+)(\</\w+\>)',
                                        '\\1\\4\\6',
                                        dom_tree.toprettyxml(encoding='utf-8')
                                        )
                                 )
                          )

        finally:
            outfile.close()

    @property
    def uses_attributes(self):
        """Property returning a Boolean value indicating that values are stored in tag attributes rather than as text
        """
        return self._uses_attributes


def main():
    # Test data from file LS7_ETM_OTH_P51_GALPGS01_092_085_20100315/scene01/LE7_20100315_092_085_L1T.xml
    TESTXML = """<?xml version="1.0" encoding="UTF-8" ?>
<EODS_DATASET>
    <MDRESOURCE>
        <MDFILEID></MDFILEID>
            <FILESIZE>842</FILESIZE>
            <RESOLUTIONINMETRES>50.000000000000 25.000000000000 12.500000000000</RESOLUTIONINMETRES>
            <CONSTRAINTID></CONSTRAINTID>
            <RESOURCESTATUS>COMPLETED</RESOURCESTATUS>
            <KEYWORDS></KEYWORDS>
            <TOPICCATEGORIES></TOPICCATEGORIES>
            <CITATION>
        <TITLE></TITLE>
        <ALTERNATETITLE>Landsat7 RCC-L1T 0920852010074</ALTERNATETITLE>
        <DATE>20111025T05:14:51</DATE>
        <DATETYPE>creation</DATETYPE>
                <EDITION></EDITION>
                <EDITIONDATE></EDITIONDATE>
                <ENTEREDBYRESPONSIBLEPARTY>
                    <INDIVIDUALNAME></INDIVIDUALNAME>
                    <ORGANISATIONNAME></ORGANISATIONNAME>
                    <POSITIONNAME></POSITIONNAME>
                    <ROLE></ROLE>
                </ENTEREDBYRESPONSIBLEPARTY>
                <UPDATEDBYRESPONSIBLEPARTY>
                    <INDIVIDUALNAME></INDIVIDUALNAME>
                    <ORGANISATIONNAME></ORGANISATIONNAME>
                    <POSITIONNAME></POSITIONNAME>
                    <ROLE></ROLE>
                </UPDATEDBYRESPONSIBLEPARTY>
                <OTHERCITATIONDETAILS></OTHERCITATIONDETAILS>
            </CITATION>
            <MDSTANDARDNAME>ANZLIC Metadata Profile: An Australian/New Zealand Profile of AS/NZS ISO 19115:2005, Geographic information - Metadata</MDSTANDARDNAME>
            <MDSTANDARDVERSION>1.1</MDSTANDARDVERSION>
            <PARENTID></PARENTID>
            <DATALANGUAGE>eng</DATALANGUAGE>
            <MDCONTACT>
        <INDIVIDUALNAME></INDIVIDUALNAME>
        <ORGANISATIONNAME></ORGANISATIONNAME>
        <POSITIONNAME></POSITIONNAME>
        <ROLE></ROLE>
            </MDCONTACT>
            <RESOURCETYPE>Processed Image</RESOURCETYPE>
            <CHARACTERSETCODE>utf8</CHARACTERSETCODE>
            <ABSTRACT>Landsat7 RCC-L1T</ABSTRACT>
            <PURPOSE></PURPOSE>
            <CREDIT></CREDIT>
            <HEIRACHYLEVEL>dataset</HEIRACHYLEVEL>
            <HEIRARCHYLEVELNAME></HEIRARCHYLEVELNAME>
            <ENVIRONMENTDESCRIPTION></ENVIRONMENTDESCRIPTION>
            <SPATIALREPRESENTATIONTYPE></SPATIALREPRESENTATIONTYPE>
            <SUPPLEMENTARYINFORMATION></SUPPLEMENTARYINFORMATION>
            <FORMAT>
        <FORMATNAME>FASTL7A</FORMATNAME>
        <FORMATVERSION></FORMATVERSION>
            </FORMAT>
        </MDRESOURCE>
    <EXEXTENT>
        <COORDINATEREFERENCESYSTEM></COORDINATEREFERENCESYSTEM>
        <EXTENTTYPE></EXTENTTYPE>
        <EXTENTDESCRIPTION></EXTENTDESCRIPTION>
        <UL_LAT>-35.0700000</UL_LAT>
        <UL_LONG>144.8700000</UL_LONG>
        <UR_LAT>-35.0700000</UR_LAT>
        <UR_LONG>147.5900000</UR_LONG>
        <LR_LAT>-37.0000000</LR_LAT>
        <LR_LONG>147.5900000</LR_LONG>
        <LL_LAT>-37.0000000</LL_LAT>
        <LL_LONG>144.8700000</LL_LONG>
        <WEST_BLONG>144.8700000</WEST_BLONG>
        <EAST_BLONG>147.5900000</EAST_BLONG>
        <NORTH_BLAT>-37.0000000</NORTH_BLAT>
        <SOUTH_BLAT>-35.0700000</SOUTH_BLAT>
        <TEMPORALEXTENTFROM>20100315 23:54:58</TEMPORALEXTENTFROM>
        <TEMPORALEXTENTTO>20100315 23:55:25</TEMPORALEXTENTTO>
        <VERTICALEXTENTMAX></VERTICALEXTENTMAX>
        <VERTICALEXTENTMIN></VERTICALEXTENTMIN>
        <VERTICALEXTENTUOM></VERTICALEXTENTUOM>
        <VERTICALEXTENTCRS></VERTICALEXTENTCRS>
        <VERTICALEXTENTDATUM></VERTICALEXTENTDATUM>
        <SCENECENTRELAT>-36.0459</SCENECENTRELAT>
        <SCENECENTRELONG>146.2550</SCENECENTRELONG>
  <TIMESERIESCOMPOSITINGINTERVAL></TIMESERIESCOMPOSITINGINTERVAL>
    </EXEXTENT>
    <DATAQUALITY>
        <SCOPELEVEL>dataset</SCOPELEVEL>
        <SCOPELEVELDESCRIPTION></SCOPELEVELDESCRIPTION>
        <LINEAGE>
            <STATEMENT>Resampling=CC,RadiometricCorrection=CPF,Orientation=NUP,LAM_version=6.2.1,LACS_version=6.3.0,LPS_version=8.2.1,LPGS_version=11.4.0</STATEMENT>
          <PROCESSINGSTEP>
   <ALGORITHMCITATION>
     <TITLE>Pinkmatter Landsat Processor</TITLE>
     <EDITION>3.2.1518</EDITION>
   </ALGORITHMCITATION>
   </PROCESSINGSTEP>
            <LINEAGESOURCE>
                <MDFILEID></MDFILEID>
          <SOURCERESOURCEID></SOURCERESOURCEID>
                <SOURCECITATION>
                    <TITLE></TITLE>
                    <DATE>20111025T05:14:51</DATE>
                    <DATETYPE>creation</DATETYPE>
               <ENTEREDBYRESPONSIBLEPARTY>
                        <INDIVIDUALNAME></INDIVIDUALNAME>
                        <ORGANISATIONNAME></ORGANISATIONNAME>
                        <POSITIONNAME></POSITIONNAME>
                        <ROLE></ROLE>
               </ENTEREDBYRESPONSIBLEPARTY>
                </SOURCECITATION>
                <DESCRIPTION></DESCRIPTION>
            <SOURCEREFERENCESYSTEM></SOURCEREFERENCESYSTEM>
            <SOURCESCALE></SOURCESCALE>
            <SOURCESTEP></SOURCESTEP>
            <PROCESSINGLEVEL></PROCESSINGLEVEL>
            </LINEAGESOURCE>
        </LINEAGE>
        <DQELEMENT>
            <MEASURENAME>Automatically Generated Report</MEASURENAME>
            <QUANTATIVEVALUE>1</QUANTATIVEVALUE>
            <QUANTATIVEVALUEUNIT>text</QUANTATIVEVALUEUNIT>
        </DQELEMENT>
    </DATAQUALITY>
    <IMAGEDESCRIPTION>
        <ILLUMINATIONELEVATIONANGLE>41.5610274</ILLUMINATIONELEVATIONANGLE>
        <ILLUMINATIONELEVATIONAZIMUTH>53.7730377</ILLUMINATIONELEVATIONAZIMUTH>
         <VIEWINGINCIDENCEANGLEXTRACK></VIEWINGINCIDENCEANGLEXTRACK>
        <VIEWINGINCIDENCEANGLELONGTRACK></VIEWINGINCIDENCEANGLELONGTRACK>
        <CLOUDCOVERPERCENTAGE></CLOUDCOVERPERCENTAGE>
        <CLOUDCOVERDETAILS></CLOUDCOVERDETAILS>
        <SENSOROPERATIONMODE>BUMPER</SENSOROPERATIONMODE>
        <BANDGAIN>HHHLHLHHL</BANDGAIN>
        <BANDSAVAILABLE>123456678</BANDSAVAILABLE>
        <SATELLITEREFERENCESYSTEM_X>092</SATELLITEREFERENCESYSTEM_X>
        <SATELLITEREFERENCESYSTEM_Y>085</SATELLITEREFERENCESYSTEM_Y>
        <IMAGECONDITION></IMAGECONDITION>
        <PROCESSINGTYPECD>L1T</PROCESSINGTYPECD>
      <BEARING></BEARING>
    </IMAGEDESCRIPTION>
    <BROWSEGRAPHIC>
        <FILENAME></FILENAME>
        <FILEDESCRIPTION></FILEDESCRIPTION>
        <FILETYPE></FILETYPE>
        <SAMPLEPIXELRESOLUTION></SAMPLEPIXELRESOLUTION>
        <BLUEBAND></BLUEBAND>
        <GREENORGREYBAND></GREENORGREYBAND>
        <REDBAND></REDBAND>
    </BROWSEGRAPHIC>
    <ACQUISITIONINFORMATION>
        <PLATFORMNAME>Landsat-7</PLATFORMNAME>
        <INSTRUMENTNAME>ETM+</INSTRUMENTNAME>
        <INSTRUMENTYPE>Multi-spectral</INSTRUMENTYPE>
        <MISSIONNAME>Landsat Data Continuity Mission (LDCM)</MISSIONNAME>
        <EVENT>
            <TIME>2010-03-15</TIME>
            <AOS>2010-03-15</AOS>
            <LOS>2010-03-15</LOS>
            <ORBITNUMBER></ORBITNUMBER>
            <CYCLENUMBER></CYCLENUMBER>
            <PASSSTATUS></PASSSTATUS>
            <NUMBERSCENESINPASS></NUMBERSCENESINPASS>
            <COLLECTIONSITE></COLLECTIONSITE>
            <ANTENNA></ANTENNA>
            <HEADING></HEADING>
            <SEQUENCE></SEQUENCE>
            <TRIGGER></TRIGGER>
            <CONTEXT></CONTEXT>
        </EVENT>
    </ACQUISITIONINFORMATION>
    <GRIDSPATIALREPRESENTATION>
        <NUMBEROFDIMENSIONS>2</NUMBEROFDIMENSIONS>
        <TRANSFORMATIONPARAMETERAVAILABILITY></TRANSFORMATIONPARAMETERAVAILABILITY>
        <CELLGEOMETRY></CELLGEOMETRY>
        <DIMENSION_X>
            <NAME>sample</NAME>
            <SIZE>5441 10881 21761</SIZE>
            <RESOLUTION>0.000500 0.000250 0.000125</RESOLUTION>
        </DIMENSION_X>
        <DIMENSION_Y>
            <NAME>line</NAME>
            <SIZE>3861 7721 15441</SIZE>
            <RESOLUTION>0.000500 0.000250 0.000125</RESOLUTION>
        </DIMENSION_Y>
        <GEORECTIFIED>
            <CHECKPOINTAVAILABILITY></CHECKPOINTAVAILABILITY>
            <CHECKPOINTDESCRIPTION></CHECKPOINTDESCRIPTION>
            <POINTINPIXEL></POINTINPIXEL>
              <GEOREFULPOINT_X></GEOREFULPOINT_X>
              <GEOREFULPOINT_Y></GEOREFULPOINT_Y>
              <GEOREFULPOINT_Z></GEOREFULPOINT_Z>
              <GEOREFURPOINT_X></GEOREFURPOINT_X>
              <GEOREFURPOINT_Y></GEOREFURPOINT_Y>
              <GEOREFURPOINT_Z></GEOREFURPOINT_Z>
              <GEOREFLLPOINT_X></GEOREFLLPOINT_X>
              <GEOREFLLPOINT_Y></GEOREFLLPOINT_Y>
              <GEOREFLLPOINT_Z></GEOREFLLPOINT_Z>
              <GEOREFLRPOINT_X></GEOREFLRPOINT_X>
              <GEOREFLRPOINT_Y></GEOREFLRPOINT_Y>
              <GEOREFLRPOINT_Z></GEOREFLRPOINT_Z>
              <CENTREPOINT_X></CENTREPOINT_X>
              <CENTREPOINT_Y></CENTREPOINT_Y>
              <ELLIPSOID>WGS84</ELLIPSOID>
              <DATUM>WGS84</DATUM>
              <ZONE></ZONE>
              <PROJECTION>EQR</PROJECTION>
            <COORDINATEREFERENCESYSTEM></COORDINATEREFERENCESYSTEM>
        </GEORECTIFIED>
    </GRIDSPATIALREPRESENTATION>
</EODS_DATASET>"""

    # Instantiate empty MTLMetadata object and parse test string (strip all EOLs first)
    xml_object = XMLMetadata()
    xml_object._populate_dict_from_node(xml.dom.minidom.parseString(TESTXML.translate(None, '\n')),
                                        xml_object.metadata_dict)
    assert xml_object.metadata_dict, 'No metadata_dict created'
    assert xml_object.tree_to_list(), 'Unable to create list from metadata_dict'
    assert xml_object.get_metadata('EODS_DATASET,ACQUISITIONINFORMATION,PLATFORMNAME'.split(',')), 'Unable to find value for key L1_METADATA_FILE,PRODUCT_METADATA,SPACECRAFT_ID'
    assert xml_object.get_metadata('...,PLATFORMNAME'.split(',')), 'Unable to find value for key ...,SPACECRAFT_ID'
    assert not xml_object.get_metadata('RUBBERCHICKEN'.split(',')), 'Found nonexistent key RUBBERCHICKEN'
    xml_object.set_metadata_node('EODS_DATASET,ACQUISITIONINFORMATION,PLATFORMNAME'.split(','), 'Rubber Chicken')
    assert xml_object.get_metadata('...,PLATFORMNAME'.split(',')), 'Unable to change ...,SPACECRAFT_ID to "Rubber Chicken"'
    xml_object.merge_metadata_dicts({'RUBBERCHICKEN': 'Rubber Chicken'}, xml_object.metadata_dict)
    assert xml_object.get_metadata('RUBBERCHICKEN'.split(',')), 'Unable to find value for key RUBBERCHICKEN'
    xml_object.delete_metadata('RUBBERCHICKEN'.split(','))
    assert not xml_object.get_metadata('RUBBERCHICKEN'.split(',')), 'Found value for key RUBBERCHICKEN'
    print xml_object.tree_to_list()
if __name__ == '__main__':
    main()
