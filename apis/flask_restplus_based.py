from flask import Flask, jsonify, request
from flask_restplus import Resource, Api, fields
from normalizer.normalizer import mutalyzer3, get_reference_models
from hgvsparser.hgvs_parser import HgvsParser
from hgvsparser import to_model

app = Flask(__name__)
api = Api(app, version='1.0', title='Mutalyzer3 API')

ns = api.namespace('/')


@ns.route("/syntax_check/<string:hgvs_description>")
class SyntaxCheck(Resource):
    def get(self, hgvs_description):
        """Takes a variant description as input and checks
        whether its syntax is correct."""
        parser = HgvsParser()
        parse_tree = parser.parse(hgvs_description)
        if parse_tree:
            return jsonify('OK')
        else:
            return jsonify('not OK')


@ns.route("/description_to_model/<string:hgvs_description>")
class DescriptionToModel(Resource):
    def get(self, hgvs_description):
        """Converts a variant description to its data model."""
        parser = HgvsParser()
        parse_tree = parser.parse(hgvs_description)
        return jsonify(to_model.convert(parse_tree))


@ns.route("/name_check/<string:hgvs_description>")
class NameCheck(Resource):
    def get(self, hgvs_description):
        """Takes a variant description as input and checks
        whether it is correct."""
        return jsonify(mutalyzer3(hgvs_description))


@ns.route("/reference_models/<string:hgvs_description>")
class ReferenceModels(Resource):
    def get(self, hgvs_description):
        """Takes a variant description as input and retrieves
        the reference models for all the references present in
        the description."""
        parser = HgvsParser()
        parse_tree = parser.parse(hgvs_description)
        description_model = to_model.convert(parse_tree)
        return jsonify(get_reference_models(description_model))


app.run(port=5001, debug=True)
