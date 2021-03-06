/*
 * Copyright 2016 EMBL - European Bioinformatics Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package uk.ac.ebi.eva.dbmigration.mongodb;

import com.mongodb.MongoClient;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;
import com.mongodb.client.MongoDatabase;
import com.mongodb.client.model.IndexOptions;

import org.bson.Document;
import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static junit.framework.TestCase.assertNull;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.ANNOT_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.CACHE_VERSION_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.CHROMOSOME_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.CONSEQUENCE_TYPE_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.DEFAULT_VERSION_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.END_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.ID_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.POLYPHEN_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.SCORE_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.SIFT_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.SO_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.START_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.VEP_VERSION_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.XREFS_FIELD;
import static uk.ac.ebi.eva.dbmigration.mongodb.ExtractAnnotationFromVariant.XREF_ID_FIELD;

/**
 * Test {@link ExtractAnnotationFromVariant}
 */
public class ExtractAnnotationFromVariantTest {

    private static final String VEP_VERSION = "88";

    private static final String CACHE_VERSION = "90";

    private static final String VARIANT_COLLECTION_NAME = "variantsCollection";

    private static final String ANNOTATION_COLLECTION_NAME = "annotationCollection";

    private static final String ANNOTATION_METADATA_COLLECTION_NAME = "annotationMetadataCollection";

    private static final String READ_PREFERENCE = "primary";

    private static final String IDS_FIELD = "ids";

    private static final String FILES_FIELD = "files";

    private final static String FILEID_FIELD = "fid";

    private final static String STUDYID_FIELD = "sid";

    private static final String ANNOTATION_FIELD = "annot";

    private static final String LEGACY_ANNOTATION_CT_SO_INDEX = "annot.ct.so";

    private static final String LEGACY_ANNOTATION_XREF_ID_INDEX = "annot.xrefs.id";

    private ExtractAnnotationFromVariant extractAnnotationFromVariant;

    private MongoClient mongoClient;

    private static final Map<String, String> databaseMap = Stream.of(new String[][]{
            {"VARIANT_WITHOUT_ANNOTATION", "variantWithoutAnnotation"},
            {"VARIANT_WITH_ANNOTATION", "variantWithAnnotation"},
            {"DB_FOR_ANNOTATION_METADATA_CHECK", "DBForAnnotationMetadataCheck"},
            {"DB_FOR_INDEXES_CHECK", "DBForIndexesCheck"},
            {"DEFAULT_ANNOTATION", "defaultAnnotation"}
    }).collect(Collectors.toMap(data -> data[0], data -> data[1]));

    @Rule
    public final ExpectedException exception = ExpectedException.none();

    @Before
    public void setUp() throws Exception {
        extractAnnotationFromVariant = new ExtractAnnotationFromVariant();
        mongoClient = new MongoClient();
    }

    @After
    public void tearDown() throws Exception {
        List<String> mongoDatabases = new ArrayList<>();
        mongoClient.listDatabaseNames().iterator().forEachRemaining(mongoDatabases::add);
        for(String usedDatabaseInTest: databaseMap.values()) {
            if(mongoDatabases.contains(usedDatabaseInTest)) {
                mongoClient.dropDatabase(usedDatabaseInTest);
            }
        }
        mongoClient.close();
    }

    @Test
    public void variantWithoutAnnotationShouldNotChange() {
        // given
        String dbName = databaseMap.get("VARIANT_WITHOUT_ANNOTATION");

        Properties properties = new Properties();
        properties.put(DatabaseParameters.VEP_VERSION, VEP_VERSION);
        properties.put(DatabaseParameters.VEP_CACHE_VERSION, CACHE_VERSION);
        properties.put(DatabaseParameters.DB_NAME, dbName);
        properties.put(DatabaseParameters.DB_COLLECTIONS_VARIANTS_NAME, VARIANT_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATIONS_NAME, ANNOTATION_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATION_METADATA_NAME, ANNOTATION_METADATA_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_READ_PREFERENCE, READ_PREFERENCE);
        DatabaseParameters databaseParameters = new DatabaseParameters();
        databaseParameters.load(properties);
        ExtractAnnotationFromVariant.setDatabaseParameters(databaseParameters);

        MongoDatabase database = mongoClient.getDatabase(dbName);
        MongoCollection<Document> variantsCollection = database.getCollection(VARIANT_COLLECTION_NAME);
        MongoCollection<Document> annotationCollection = database.getCollection(ANNOTATION_COLLECTION_NAME);


        Document variantWithoutAnnot = Document.parse(VariantData.VARIANT_WITHOUT_ANNOT);
        variantsCollection.insertOne(variantWithoutAnnot);

        createLegacyIndexes(variantsCollection);

        // when
        extractAnnotationFromVariant.migrateAnnotation(database);
        extractAnnotationFromVariant.reduceAnnotationFromVariants(database);

        // then
        try (MongoCursor<Document> variantCursor = variantsCollection.find().iterator()) {
            while (variantCursor.hasNext()) {
                Document variantObj = variantCursor.next();
                assertNull(variantObj.get(ANNOT_FIELD));
            }
        }

        assertEquals(0, variantsCollection.count(new Document(ANNOT_FIELD, new Document("$exists", true))));
        assertEquals(0, annotationCollection.count());
    }

    @Test
    public void variantWithAnnotationShouldMigrate() {
        // given
        String dbName = databaseMap.get("VARIANT_WITH_ANNOTATION");

        Properties properties = new Properties();
        properties.put(DatabaseParameters.VEP_VERSION, VEP_VERSION);
        properties.put(DatabaseParameters.VEP_CACHE_VERSION, CACHE_VERSION);
        properties.put(DatabaseParameters.DB_NAME, dbName);
        properties.put(DatabaseParameters.DB_COLLECTIONS_VARIANTS_NAME, VARIANT_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATIONS_NAME, ANNOTATION_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATION_METADATA_NAME, ANNOTATION_METADATA_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_READ_PREFERENCE, READ_PREFERENCE);
        DatabaseParameters databaseParameters = new DatabaseParameters();
        databaseParameters.load(properties);
        ExtractAnnotationFromVariant.setDatabaseParameters(databaseParameters);

        MongoDatabase database = mongoClient.getDatabase(dbName);
        MongoCollection<Document> variantsCollection = database.getCollection(VARIANT_COLLECTION_NAME);
        MongoCollection<Document> annotationCollection = database.getCollection(ANNOTATION_COLLECTION_NAME);

        Document variantWithAnnot = Document.parse(VariantData.VARIANT_WITH_ANNOT_1);
        variantsCollection.insertOne(variantWithAnnot);

        Document originalVariant = variantsCollection.find().first();
        Document originalAnnotField = (Document)originalVariant.get(ANNOT_FIELD);

        createLegacyIndexes(variantsCollection);

        // when
        extractAnnotationFromVariant.migrateAnnotation(database);

        // then
        assertEquals(1, annotationCollection.count());

        Document annotation = annotationCollection.find().first();

        String versionSuffix = "_" + databaseParameters.getVepVersion() + "_" + databaseParameters.getVepCacheVersion();
        assertEquals(originalVariant.get(ID_FIELD) + versionSuffix, annotation.get(ID_FIELD));
        assertEquals(originalVariant.get(CHROMOSOME_FIELD), annotation.get(CHROMOSOME_FIELD));
        assertEquals(originalVariant.get(START_FIELD), annotation.get(START_FIELD));
        assertEquals(originalVariant.get(END_FIELD), annotation.get(END_FIELD));
        assertEquals(VEP_VERSION, annotation.get(VEP_VERSION_FIELD));
        assertEquals(CACHE_VERSION, annotation.get(CACHE_VERSION_FIELD));
        assertEquals(originalAnnotField.get(CONSEQUENCE_TYPE_FIELD), annotation.get(CONSEQUENCE_TYPE_FIELD));
        assertEquals(originalAnnotField.get(XREFS_FIELD), annotation.get(XREFS_FIELD));
    }

    @Test
    public void variantWithAnnotationShouldKeepSomeFields() {
        // given
        String dbName = databaseMap.get("VARIANT_WITH_ANNOTATION");

        Properties properties = new Properties();
        properties.put(DatabaseParameters.VEP_VERSION, VEP_VERSION);
        properties.put(DatabaseParameters.VEP_CACHE_VERSION, CACHE_VERSION);
        properties.put(DatabaseParameters.DB_NAME, dbName);
        properties.put(DatabaseParameters.DB_COLLECTIONS_VARIANTS_NAME, VARIANT_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATIONS_NAME, ANNOTATION_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATION_METADATA_NAME, ANNOTATION_METADATA_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_READ_PREFERENCE, READ_PREFERENCE);
        DatabaseParameters databaseParameters = new DatabaseParameters();
        databaseParameters.load(properties);
        ExtractAnnotationFromVariant.setDatabaseParameters(databaseParameters);

        MongoDatabase database = mongoClient.getDatabase(dbName);
        MongoCollection<Document> variantsCollection = database.getCollection(VARIANT_COLLECTION_NAME);

        Document variantWithAnnot = Document.parse(VariantData.VARIANT_WITH_ANNOT_2);
        variantsCollection.insertOne(variantWithAnnot);

        Document originalVariant = variantsCollection.find().first();
        Document originalAnnotField = (Document) originalVariant.get(ANNOT_FIELD);

        createLegacyIndexes(variantsCollection);

        // when
        extractAnnotationFromVariant.reduceAnnotationFromVariants(database);

        // then
        assertEquals(1, variantsCollection.count());

        Document variant = variantsCollection.find().first();
        List newAnnotField = (List) variant.get(ANNOT_FIELD);
        Document newAnnotElement = (Document) newAnnotField.get(0);

        assertEquals(VEP_VERSION, newAnnotElement.get(VEP_VERSION_FIELD));
        assertEquals(CACHE_VERSION, newAnnotElement.get(CACHE_VERSION_FIELD));

        ArrayList<Integer> so = (ArrayList<Integer>) newAnnotElement.get(SO_FIELD);
        Set<Integer> expectedSo = computeSo(originalAnnotField);
        assertEquals(expectedSo.size(), so.size());
        assertEquals(expectedSo, new TreeSet<>(so));

        ArrayList<String> xrefs = (ArrayList<String>) newAnnotElement.get(XREFS_FIELD);
        Set<String> expectedXref = computeXref(originalAnnotField);
        assertEquals(expectedXref.size(), xrefs.size());
        assertEquals(expectedXref, new TreeSet<>(xrefs));

        assertEquals(computeSift(originalAnnotField), newAnnotElement.get(SIFT_FIELD));

        assertEquals(computePolyphen(originalAnnotField), newAnnotElement.get(POLYPHEN_FIELD));
    }

    private Set<Integer> computeSo(Document originalAnnotField) {
        Set<Integer> soSet = new TreeSet<>();

        List<Document> cts = (List<Document>) originalAnnotField.get(CONSEQUENCE_TYPE_FIELD);
        for (Document ct : cts) {
            soSet.addAll(((List<Integer>) ct.get(SO_FIELD)));
        }

        return soSet;
    }

    private Set<String> computeXref(Document originalAnnotField) {
        Set<String> xrefSet = new TreeSet<>();

        List<Document> cts = (List<Document>) originalAnnotField.get(XREFS_FIELD);
        for (Document ct : cts) {
            xrefSet.add(((String) ct.get(XREF_ID_FIELD)));
        }

        return xrefSet;
    }

    private List<Double> computeSift(Document originalAnnotField) {
        Double min = Double.POSITIVE_INFINITY;
        Double max = Double.NEGATIVE_INFINITY;

        List<Document> cts = (List<Document>) originalAnnotField.get(CONSEQUENCE_TYPE_FIELD);
        for (Document ct : cts) {
            Document document = (Document) ct.get(SIFT_FIELD);
            if (document != null) {
                Double score = (Double) document.get(SCORE_FIELD);

                min = Math.min(min, score);
                max = Math.max(max, score);
            }
        }

        return Arrays.asList(min, max);
    }

    private List<Double> computePolyphen(Document originalAnnotField) {
        Double min = Double.POSITIVE_INFINITY;
        Double max = Double.NEGATIVE_INFINITY;

        List<Document> cts = (List<Document>) originalAnnotField.get(CONSEQUENCE_TYPE_FIELD);
        for (Document ct : cts) {
            Document document = (Document) ct.get(POLYPHEN_FIELD);
            if (document != null) {
                Double score = (Double) document.get(SCORE_FIELD);

                min = Math.min(min, score);
                max = Math.max(max, score);
            }
        }
        return Arrays.asList(min, max);
    }

    @Test
    public void metadataShouldBeUpdated() throws Exception {
        // given
        String dbName = databaseMap.get("DB_FOR_ANNOTATION_METADATA_CHECK");

        Properties properties = new Properties();
        properties.put(DatabaseParameters.VEP_VERSION, VEP_VERSION);
        properties.put(DatabaseParameters.VEP_CACHE_VERSION, CACHE_VERSION);
        properties.put(DatabaseParameters.DB_NAME, dbName);
        properties.put(DatabaseParameters.DB_COLLECTIONS_VARIANTS_NAME, VARIANT_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATIONS_NAME, ANNOTATION_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATION_METADATA_NAME, ANNOTATION_METADATA_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_READ_PREFERENCE, READ_PREFERENCE);
        DatabaseParameters databaseParameters = new DatabaseParameters();
        databaseParameters.load(properties);
        ExtractAnnotationFromVariant.setDatabaseParameters(databaseParameters);

        MongoDatabase database = mongoClient.getDatabase(dbName);
        MongoCollection<Document> variantsCollection = database.getCollection(VARIANT_COLLECTION_NAME);
        MongoCollection<Document> annotationMetadataCollection = database.getCollection(
                ANNOTATION_METADATA_COLLECTION_NAME);

        Document variantWithAnnot = Document.parse(VariantData.VARIANT_WITH_ANNOT_2);
        variantsCollection.insertOne(variantWithAnnot);

        // when
        extractAnnotationFromVariant.updateAnnotationMetadata(database);

        // then
        assertEquals(1, annotationMetadataCollection.count());
        assertEquals(VEP_VERSION + "_" + CACHE_VERSION, annotationMetadataCollection.find().first().get(ID_FIELD));
        assertEquals(VEP_VERSION, annotationMetadataCollection.find().first().get(VEP_VERSION_FIELD));
        assertEquals(CACHE_VERSION, annotationMetadataCollection.find().first().get(CACHE_VERSION_FIELD));
    }

    @Test
    public void indexesShouldBeCreated() throws Exception {
        // given
        String dbName = databaseMap.get("DB_FOR_INDEXES_CHECK");

        Properties properties = new Properties();
        properties.put(DatabaseParameters.VEP_VERSION, VEP_VERSION);
        properties.put(DatabaseParameters.VEP_CACHE_VERSION, CACHE_VERSION);
        properties.put(DatabaseParameters.DB_NAME, dbName);
        properties.put(DatabaseParameters.DB_COLLECTIONS_VARIANTS_NAME, VARIANT_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATIONS_NAME, ANNOTATION_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATION_METADATA_NAME, ANNOTATION_METADATA_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_READ_PREFERENCE, READ_PREFERENCE);
        DatabaseParameters databaseParameters = new DatabaseParameters();
        databaseParameters.load(properties);
        ExtractAnnotationFromVariant.setDatabaseParameters(databaseParameters);

        MongoDatabase database = mongoClient.getDatabase(dbName);
        MongoCollection<Document> variantsCollection = database.getCollection(VARIANT_COLLECTION_NAME);
        MongoCollection<Document> annotationsCollection = database.getCollection(ANNOTATION_COLLECTION_NAME);

        Document variantWithAnnot = Document.parse(VariantData.VARIANT_WITH_ANNOT_2);
        variantsCollection.insertOne(variantWithAnnot);

        createLegacyIndexes(variantsCollection);

        // when
        extractAnnotationFromVariant.migrateAnnotation(database);
        extractAnnotationFromVariant.dropIndexes(database);
        extractAnnotationFromVariant.reduceAnnotationFromVariants(database);
        extractAnnotationFromVariant.updateAnnotationMetadata(database);
        extractAnnotationFromVariant.createIndexes(database);

        // then
        ArrayList<Document> variantsIndexes = variantsCollection.listIndexes().into(new ArrayList<>());
        assertEquals(6, variantsIndexes.size());
        assertIndexNameExists(variantsIndexes, ANNOT_FIELD + "." + SO_FIELD + "_1");
        assertIndexNameExists(variantsIndexes, ANNOT_FIELD + "." + XREFS_FIELD + "_1");

        ArrayList<Document> annotationsIndexes = annotationsCollection.listIndexes().into(new ArrayList<>());
        assertEquals(4, annotationsIndexes.size());
        assertIndexNameExists(annotationsIndexes, CONSEQUENCE_TYPE_FIELD + "." + SO_FIELD + "_1");
        assertIndexNameExists(annotationsIndexes, XREFS_FIELD + "." + XREF_ID_FIELD + "_1");
        assertIndexNameExists(annotationsIndexes,
                              String.join("_", CHROMOSOME_FIELD, "1", START_FIELD, "1", END_FIELD, "1"));
    }

    private void createLegacyIndexes(MongoCollection<Document> variantsCollection) {
        IndexOptions background = new IndexOptions().background(true);
        variantsCollection.createIndex(
                new Document(CHROMOSOME_FIELD, 1).append(START_FIELD, 1).append(END_FIELD, 1),
                background);

        variantsCollection.createIndex(new Document(IDS_FIELD, 1), background);

        String filesStudyIdField = String.format("%s.%s", FILES_FIELD, STUDYID_FIELD);
        String filesFileIdField = String.format("%s.%s", FILES_FIELD, FILEID_FIELD);
        variantsCollection.createIndex(new Document(filesStudyIdField, 1).append(filesFileIdField, 1), background);

        variantsCollection.createIndex(new Document(LEGACY_ANNOTATION_CT_SO_INDEX, 1), background);
        variantsCollection.createIndex(new Document(LEGACY_ANNOTATION_XREF_ID_INDEX, 1), background);
    }

    private void assertIndexNameExists(ArrayList<Document> variantsIndexes, String indexName) {
        boolean indexFound = false;
        for (Document variantsIndex : variantsIndexes) {
            String name = variantsIndex.get("name", String.class);
            if (name != null) {
                if (name.equals(indexName)) {
                    indexFound = true;
                }
            }
        }
        assertTrue(indexFound);
    }

    @Test
    public void testDefaultAnnotationVersion() throws Exception {
        String dbName = databaseMap.get("DEFAULT_ANNOTATION");

        Properties properties = new Properties();
        properties.put(DatabaseParameters.VEP_VERSION, VEP_VERSION);
        properties.put(DatabaseParameters.VEP_CACHE_VERSION, CACHE_VERSION);
        properties.put(DatabaseParameters.DB_NAME, dbName);
        properties.put(DatabaseParameters.DB_COLLECTIONS_VARIANTS_NAME, VARIANT_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATIONS_NAME, ANNOTATION_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATION_METADATA_NAME, ANNOTATION_METADATA_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_READ_PREFERENCE, READ_PREFERENCE);
        DatabaseParameters databaseParameters = new DatabaseParameters();
        databaseParameters.load(properties);
        ExtractAnnotationFromVariant.setDatabaseParameters(databaseParameters);

        MongoDatabase mongoDatabase = mongoClient.getDatabase(dbName);
        MongoCollection<Document> collection = mongoDatabase.getCollection(ANNOTATION_METADATA_COLLECTION_NAME);

        collection.insertOne(buildAnnotationMetadataDocument("79", "78"));
        collection.insertOne(buildAnnotationMetadataDocument("80", "82"));
        collection.insertOne(buildAnnotationMetadataDocument(VEP_VERSION, CACHE_VERSION));
        collection.insertOne(buildAnnotationMetadataDocument("98", "99"));

        extractAnnotationFromVariant.addDefaultVersion(mongoDatabase);

        for (Document document : collection.find()) {
            if (document.get(VEP_VERSION_FIELD).equals(VEP_VERSION)) {
                assertTrue((Boolean) document.get(DEFAULT_VERSION_FIELD));
            } else {
                assertFalse((Boolean) document.get(DEFAULT_VERSION_FIELD));
            }
        }
    }

    @Test
    public void testDefaultAnnotationVersionNotFound() throws Exception {
        String dbName = databaseMap.get("DEFAULT_ANNOTATION");

        Properties properties = new Properties();
        properties.put(DatabaseParameters.VEP_VERSION, VEP_VERSION);
        properties.put(DatabaseParameters.VEP_CACHE_VERSION, CACHE_VERSION);
        properties.put(DatabaseParameters.DB_NAME, dbName);
        properties.put(DatabaseParameters.DB_COLLECTIONS_VARIANTS_NAME, VARIANT_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATIONS_NAME, ANNOTATION_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_COLLECTIONS_ANNOTATION_METADATA_NAME, ANNOTATION_METADATA_COLLECTION_NAME);
        properties.put(DatabaseParameters.DB_READ_PREFERENCE, READ_PREFERENCE);
        DatabaseParameters databaseParameters = new DatabaseParameters();
        databaseParameters.load(properties);
        ExtractAnnotationFromVariant.setDatabaseParameters(databaseParameters);

        MongoDatabase mongoDatabase = mongoClient.getDatabase(dbName);
        MongoCollection<Document> collection = mongoDatabase.getCollection(ANNOTATION_METADATA_COLLECTION_NAME);

        collection.insertOne(buildAnnotationMetadataDocument("79", "78"));
        collection.insertOne(buildAnnotationMetadataDocument("80", "82"));
        collection.insertOne(buildAnnotationMetadataDocument("98", "99"));

        exception.expect(IllegalStateException.class);
        extractAnnotationFromVariant.addDefaultVersion(mongoDatabase);
    }

    private Document buildAnnotationMetadataDocument(String vepVersion, String vepCacheVersion) {
        return new Document(ID_FIELD, vepVersion + "_" + vepCacheVersion)
                .append(VEP_VERSION_FIELD, vepVersion)
                .append(CACHE_VERSION_FIELD, vepCacheVersion);
    }
}
