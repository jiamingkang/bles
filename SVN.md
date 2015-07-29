

&lt;B&gt;

1. Release New Version 

&lt;/B&gt;



-- Usage:
svn copy [-r ${revision-number}] ${SOURCE\_URL} ${DESTINATION\_URL} [-m "message describing the details of the version"}

-- Example:
svn copy [-r [r37](https://code.google.com/p/bles/source/detail?r=37)] \
https://bles.googlecode.com/svn/trunk \
https://bles.googlecode.com/svn/tags/OOD-v0.2 \
-m "OOD running version v0.2"



&lt;B&gt;


2. Check out the source code


&lt;/B&gt;



-- Usage:
svn co ${SOURCE\_URL} ${DESTINATION\_LOCAL\_PATH}

-- Example:
svn co http://bles.googlecode.com/svn/trunk topopt



&lt;B&gt;


3. Check in the new code into repository


&lt;/B&gt;



-- Usage:
svn ci [-m "messages describing changes in this version"]

-- Example:
svn ci -m "re-write algorithm"



&lt;B&gt;


4. Update the local source folder


&lt;/B&gt;



-- Usage:
svn up

-- Example:
svn up

-- Comment:
If you need a particular version, please use a revision option
(e.g. svn up -r [r36](https://code.google.com/p/bles/source/detail?r=36)). Revision number should start a character 'r'.



&lt;B&gt;


5. Add new items in the repository


&lt;/B&gt;



-- Usage:
svn add ${ITEM\_TO\_BE\_ADDED}

-- Example:
svn add CHakBetterCode.cpp

-- Comments:
'svn add' does not mean the commitment of the new item into the server. To complete the commitment, it should be followed by 'svn ci'. Please see below

svn add CHakBetterCode.cpp

svn ci -m "added the CHakBetterCode.cpp"

After that, the addition of the code is completed in the server repository.



&lt;B&gt;


6. Delete items in the repository


&lt;/B&gt;



-- Usage:
svn rm ${ITEM\_TO\_BE\_DELETED}

-- Example:
svn rm CHakLegacyCode.cpp

-- Comments:
'svn rm' does not mean the deletion of the item from the server. To complete the deletion, it should be followed by 'svn ci'. Please see below

svn rm CHakLegacyCode.cpp

svn ci -m "delete CHakLegacyCode.cpp"

After that, the deletion of the code is completed in the server repository.