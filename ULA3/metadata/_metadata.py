#!/usr/bin/env python

"""Metadata module

Author: Alex Ip (alex.ip@ga.gov.au)
"""
import pickle, logging, os

logger = logging.getLogger('root.' + __name__)

class MetadataException(Exception):
    pass

class Metadata(object):
    """Superclass of all metadata types
    Manages master dict containing all metadata trees
    """
    # Class variable holding metadata type string (e.g. 'XML', 'MTL', 'REPORT', 'TIF', 'FST')
    _metadata_type_id = None # Not set for the master metadata class
    _filename_pattern = '.*\.dat' # Default RegEx for finding metadata file.

    def __init__(self, source = None):
        """Instantiates Metadata object
        """
        self._metadata_dict = {}
        self._filename=None

        if source:
            if type(source) == dict:
                self._metadata_dict = source
            elif type(source) == str:
                self.read_file(source)

    def get_metadata(self, key_path_list=[], subtree=None):
        """Function to return the sub-dict or value in the metadata nested dict
        from a list of keys drilling down through the tree structure. Key path
        can also contain ellipsis ('...') to skip to the first found instance
        of the next key
        Returns:
            subtree dict, metadata value or None
        Side effect: Will pop values from the start of key_path_list until key is found
        """
        def find_first_key(search_key, search_dict):
            """Recursive helper function to find the first value or sub-dict for the specified
            search key when an ellipsis is used in a key path.
            """
            logger.debug('  find_first_key(%s, %s) called', repr(search_key), repr(search_dict))
            if type(search_dict) != dict:
                return None

            for key in search_dict.keys():
                if key == search_key:
                    return search_dict[key]
                else:
                    found_item = find_first_key(search_key, search_dict[key])
                    if found_item:
                        return found_item

        logger.debug('get_metadata(%s, %s) called', repr(key_path_list), repr(subtree))

        subtree = subtree or self._metadata_dict
#        assert subtree, 'Subtree must be specified'

        # Convert comma-delimited string to list if necessary
        if type(key_path_list) == str:
            key_path_list = key_path_list.split(',')

        key_path_list = list(key_path_list) # Do not modify original list (is this necessary?)
        while subtree and key_path_list:
            key = key_path_list.pop(0)
            if key == '...': # Ellipsis means skip to next key
                while key_path_list and key == '...': # Skip to next non-ellipsis key
                    key = key_path_list.pop(0) # Skip to next key
                if key == '...': # Bad input - ends in ellipsis
                    return None
                subtree = self.get_metadata(key_path_list, find_first_key(key, subtree))
            elif key:
                subtree = subtree.get(key)
                logger.debug('key = %s, value = %s', key, subtree)

        return subtree


    def delete_metadata(self, key_path_list, subtree=None):
        logger.debug('delete_metadata(%s, %s) called', repr(key_path_list), repr(subtree))
        assert key_path_list, "Key path list must be non-empty"
        _key_path_list = list(key_path_list) # Copy list to avoid side effects
        key = _key_path_list.pop()
        subtree = self.get_metadata(_key_path_list, subtree)
        assert subtree and key in subtree.keys(), repr(key_path_list) + " not found"
        del subtree[key]
        logger.debug('%s deleted', repr(key_path_list))


    def tree_to_tuples(self, subtree=None, node_name=''):
        """Recursive function to return all leaf node (key, value) pairs as a flat (un-sorted) list of tuples
        Arguments:
            subtree: nested dict to contain nodes. Defaults to full internal metadata dict
            node_name: comma-separated node path to pre-pend to child node names
        Returns:
            flat list of (<node path>, <value>) tuples
        """
        logger.debug('tree_to_tuples(%s, %s) called', repr(subtree), repr(node_name))

        subtree = subtree or self._metadata_dict
#        assert subtree, 'Subtree must be specified'

        result_list = []
        while subtree:
            subtree = subtree.copy() # Do not modify original top-level dict
            key, value = subtree.popitem()
            if node_name:
                key = node_name + ',' + key
            if type(value) == dict: # not a leaf node
                result_list += self.tree_to_tuples(value, key)
            else: # Leaf node - add key=value string to list
                result_list.append((str(key), str(value)))

        return result_list

    def tree_to_list(self, subtree=None, node_name=''):
        """Recursive function to return all leaf node (key, value) pairs as a flat (un-sorted) list of strings
        Arguments:
            subtree: nested dict to contain nodes. Defaults to full internal metadata dict
            node_name: comma-separated node path to pre-pend to child node names
        Returns:
            flat list of <node path>=<value> strings
        """
        logger.debug('tree_to_list(%s, %s) called', repr(subtree), repr(node_name))

        subtree = subtree or self._metadata_dict
#        assert subtree, 'Subtree must be specified'

        return [name + '=' + value for name, value in self.tree_to_tuples(subtree)]

    def read_file(self, filename=None):
        """Abstract function to parse a metadata file and store the results in self._metadata_dict
        Needs to be implemented for the relevant file format in all descendant classes
        Argument:
            filename: Name of metadata file to be parsed and stored. Defaults to instance value
        Returns:
            Nested dict containing metadata
        Side effect:
            Changes value of self._filename to match specified name if load succeeds
        """
        filename = filename or self._filename
        assert filename, 'Filename must be specified'

        infile = open(filename, 'rb')
        self._metadata_dict = pickle.load(infile)
        infile.close()

        self._filename = filename
        return self._metadata_dict

    def write_file(self, filename=None):
        """Abstract function write the metadata contained in self._metadata_dict to a
        file in the appropriate format.
        Needs to be implemented for the relevant file format in all descendant classes
        Argument:
            filename: Metadata file to be written
        """
        filename = filename or self._filename
        assert filename, 'Filename must be specified'

        outfile = open(filename, 'wb')
        pickle.dump(self._metadata_dict, outfile)
        outfile.close()

    def set_root_metadata_from_object(self, metadata_object):
        """Function to add the metadata belonging to another metadata object to the internal metadata dict
        Argument: Metadata object
        """
        # ToDo: Implement type checking to ensure that metadata_object is a Metadata instance
        self._metadata_dict[metadata_object.metadata_type_id] = metadata_object.metadata_dict
        return self._metadata_dict

    def merge_root_metadata_from_object(self, metadata_object, overwrite=True):
        """Function to merge the metadata belonging to another metadata object to the internal metadata dict
        Argument: Metadata object
        """
        # ToDo: Implement type checking to ensure that metadata_object is a Metadata instance
        self.merge_root_metadata(metadata_object.metadata_type_id, metadata_object.metadata_dict, overwrite)
        return self._metadata_dict

    def set_root_metadata(self, metadata, root_key=None):
        """Function to add or replace a nested dict under the specified root key in the internal metadata dict
        Arguments:
            root_key: metadata type string (e.g. 'XML', 'MTL', 'REPORT', 'TIF', 'FST')
            metadata: nested dict containing metadata tree to be added. Could also be a scalar value
        """
        if root_key:
            self._metadata_dict[root_key] = metadata
        else:
            self._metadata_dict = metadata

        return self._metadata_dict

    def merge_root_metadata(self, root_key, metadata, overwrite=True):
        """Function to merge a nested dict under the specified root key in the internal metadata dict.
        N.B: Will always overwrite existing values
        Arguments:
            root_key: metadata type string (e.g. 'XML', 'MTL', 'REPORT', 'TIF', 'FST')
            metadata: nested dict containing metadata tree to be added.
        """
        destination_tree = self._metadata_dict.get(root_key)
        if not destination_tree:
            destination_tree = {}
            self._metadata_dict[root_key] = destination_tree

        self.merge_metadata_dicts(metadata, destination_tree, overwrite)
        return self._metadata_dict

    def merge_metadata_node(self, key_path_list, metadata, overwrite=True):
        """Function to merge a nested dict at a node specified by a list of
        keys drilling down through the tree structure
        Arguments:
            key_path_list: List of keys defining the path to the node.
            metadata: Value or nested dict to graft into _metadata_dict
            overwrite: Boolean flag to enable overwriting of existing values
        """
        logger.debug('merge_metadata_node(%s, %s, %s) called', repr(key_path_list), repr(metadata), repr(overwrite))

        # Convert comma-delimited string to list if necessary
        if type(key_path_list) == str:
            key_path_list = key_path_list.split(',')

        destination_tree = self.get_metadata(key_path_list)
        assert destination_tree, 'Destination subtree dict not found'
        assert type(destination_tree) == dict, 'Destination is not a dict'

        assert '...' not in key_path_list, 'Key path must be specified explicitly (no ellipses allowed)'

        self.merge_metadata_dicts(metadata, destination_tree, True)
        return self._metadata_dict

    def set_metadata_node(self, key_path_list, metadata, overwrite=True):
        """Function to set a value or sub-dict in the metadata nested dict at a node specified by a list of
        keys drilling down through the tree structure.
        Arguments:
            key_path_list: List of keys defining the path to the node.
            metadata: Value or nested dict to graft into _metadata_dict
            overwrite: Boolean flag to enable overwriting of existing values
        N.B: Key path may NOT contain ellipses ('...')
        """
        logger.debug('set_metadata_node(%s, %s, %s) called', repr(key_path_list), repr(metadata), repr(overwrite))

        # Convert comma-delimited string to list if necessary
        if type(key_path_list) == str:
            key_path_list = key_path_list.split(',')

        assert key_path_list, 'Key path must be specified with a non-empty list'

        # Convert comma-delimited string to list if necessary
        if type(key_path_list) == str:
            key_path_list = key_path_list.split(',')

        assert '...' not in key_path_list, 'Key path must be specified explicitly (no ellipses allowed)'
        subtree = self._metadata_dict
        key_path_list = list(key_path_list) # Do not modify original list
        while type(subtree) == dict and key_path_list:
            key = key_path_list.pop(0)
            logger.debug('key = %s', key)
            if key:
                if not key_path_list: # No more levels to descend
                    if key in subtree.keys() and subtree[key] and not overwrite: # Metadata for key already exists in subtree
                        raise Exception('Unable to overwrite key ' + key)
                    else:
#                        logger.debug('  Setting subtree[key] = %s', metadata)
                        subtree[key] = metadata # Overwrite previous node
                else: # still more levels to descend
                    if key in subtree.keys(): # Existing node found (dict or value)
                        if type(subtree[key]) != dict and subtree[key] and not overwrite:
                            raise Exception('Unable to overwrite subtree ' + key)
                        else:
#                            logger.debug('  Setting subtree = %s', subtree.get(key))
                            subtree = subtree.get(key) # Descend to next level
                    else: # Key doesn't exist in subtree
                        logger.debug('  Setting subtree[key] = {}')
                        subtree[key] = {} # Create new subtree
                        subtree = subtree[key]

    def merge_metadata_dicts(self, source_tree, destination_tree, overwrite=True,
                             add_new_nodes=True, keep_existing_data=False):
        """Recursive function to copy a nested dict into another nested dict
        Arguments:
            source_tree: nested dict representing tree to copy into destination tree
            destination_tree: nested dict representing destination tree
            overwrite: Boolean flag to enable overwriting of existing values. Will raise exception on conflict
            add_new_nodes: Boolean flag to enable creation of new values
            keep_existing_data: Boolean flag to keep any non-empty data in the destination tree
        """
        if source_tree is None:
            return

        assert type(source_tree) == dict, 'Source tree must be a dict'

        for key in source_tree.keys():
            source_metadata = source_tree[key]
            dest_metadata = destination_tree.get(key)
            if type(source_metadata) == dict: # Source metadata is not a leaf node
                if dest_metadata is None: # Key doesn't exist in destination - create sub-dict
                    if not add_new_nodes:
                        logger.debug('Unable to create new node %s', key)
                        continue
                    dest_metadata = {}
                    destination_tree[key] = dest_metadata
                elif type(dest_metadata) != dict: # Destination metadata is a leaf node
                    # Overwrite leaf node with new sub-dict if possible
                    if not overwrite:
                        logger.debug('Unable to overwrite existing leaf node %s', key)
                        continue
                    dest_metadata = {}
                    destination_tree[key] = dest_metadata
                # Recursive call to copy source subtree to destination subtree
                self.merge_metadata_dicts(source_metadata, dest_metadata, overwrite,
                                          add_new_nodes, keep_existing_data)
            else: # source metadata is a leaf node
                if dest_metadata and not overwrite:
                    logger.debug('Unable to overwrite existing destination metadata for %s', key)
                    continue
                elif dest_metadata is None and not add_new_nodes:
                    logger.debug('Unable to create new node %s', key)
                    continue

                if not dest_metadata or not keep_existing_data:
                    destination_tree[key] = source_metadata
                else:
                    logger.debug('Kept existing value of node %s', key)

    @property
    def filename(self):
        """Returns filename
        """
        if self._filename:
            return os.path.abspath(self._filename)
        else:
            return None

    @property
    def metadata_type_id(self):
        """Returns metadata type ID string
        """
        return self.__class__._metadata_type_id

    @property
    def metadata_dict(self):
        """Returns metadata dict containing the full metadata tree
        """
        return self._metadata_dict

if __name__ == '__main__':
    # Need a better test here
    from ULA3.dataset import SceneDataset
    l = SceneDataset()
    l.Open('/home/alex/nbar/test_data/LS7_ETM_OTH_P51_GALPGS01_092_085_20100315')
    print l.GetMetadata_List()


